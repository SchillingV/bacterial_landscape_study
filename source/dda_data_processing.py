import os
import sys
import shutil
import pandas as pd
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import tempfile
import zipfile

# ---- This is an Mass spectrometry processing pipeline for dda data wrapped in Python for ease of use on HPC clusters and large amounts of different samples 
# ---- CONFIG: Please adjust these to your install layout when you try and run this pipeline --
# ---- The tools used here can be acquired from the GitHub repos from the respective labs --
PHILOSOPHER = (Path(os.getcwd()) / "fragpipe-23.0/tools/Philosopher/philosopher-v5.1.1").resolve()  # executable
MSFRAGGER_JAR = (Path(os.getcwd()) / "fragpipe-23.0/tools/MSFragger-4.2/MSFragger-4.2.jar").resolve()                   # jar
IONQUANT_JAR = (Path(os.getcwd()) / "fragpipe-23.0/tools/IonQuant-1.11.9/IonQuant-1.11.9.jar").resolve()                   # jar
BASE_FRAGGER_PARAMS = (Path(os.getcwd()) / "fragger.params").resolve()                              # template for all runs 
DECOY_PREFIX = "rev_"
JAVA_MEM_GB = 128  # per job 

def run(cmd, cwd=None):
    """Run single command"""
    print(f"[RUN] {' '.join(map(str, cmd))} (cwd={cwd})")
    subprocess.run(cmd, check=True, cwd=cwd)

def make_target_decoy_fasta(target_fasta: Path, out_fasta: Path, decoy_prefix: str = "rev_", contam_fasta: Path | None = None):
    """
    target+decoy FASTA to "out_fasta". Decoys reversed sequences
    with IDs prefixed by decoy_prefix
    """
    def write_record(fh, header, seq):
        fh.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + "\n")

    n_target = n_decoy = 0
    with out_fasta.open("w") as out:
        # targets + decoys
        with target_fasta.open() as fh:
            header, seq = None, []
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if header is not None and seq:
                        s = "".join(seq).replace("*","")
                        write_record(out, header, s)
                        write_record(out, f"{decoy_prefix}{header}", s[::-1])
                        n_target += 1; n_decoy += 1
                    header = line[1:].strip()
                    seq = []
                else:
                    seq.append(line.strip())
            if header is not None and seq:
                s = "".join(seq).replace("*","")
                write_record(out, header, s)
                write_record(out, f"{decoy_prefix}{header}", s[::-1])
                n_target += 1; n_decoy += 1

        # optional contaminants
        if contam_fasta and contam_fasta.exists():
            with contam_fasta.open() as fh:
                out.write(fh.read())

    print(f"[DB] Wrote target+decoy FASTA: {out_fasta} (targets={n_target}, decoys={n_decoy})")

def write_params_with_fasta(base_params: Path, out_params: Path, fasta_path: Path, num_threads: int = None):
    """
    fragger.params is cloned and adjusted 
    fasta is the database here 
    """
    lines = base_params.read_text().splitlines()
    have_db = False
    have_threads = False
    out = []
    for ln in lines:
        key = ln.split("=", 1)[0].strip() if "=" in ln else None
        if key == "database_name":
            out.append(f"database_name = {fasta_path}")
            have_db = True
        elif num_threads is not None and key == "num_threads":
            out.append(f"num_threads = {num_threads}")
            have_threads = True
        else:
            out.append(ln)

    if not have_db:
        out.append(f"database_name = {fasta_path}")
    if num_threads is not None and not have_threads:
        out.append(f"num_threads = {num_threads}")

    out_params.write_text("\n".join(out) + "\n")

def collect_msfragger_outputs(
    sample_dir: Path,
    d_folder: Path,
    side_files: str = "move",  # "move" | "delete" | "skip"
) -> tuple[str | None, str | None]:

    d_parent = d_folder.parent

    # --- 1) find & move pepXML ---
    pepxmls = sorted(d_parent.glob(f"{d_folder.stem}*.pepXML"))
    if not pepxmls:
        # very rare: try inside the .d itself
        pepxmls = sorted(d_folder.glob("*.pepXML"))
    if not pepxmls:
        print(f"[WARN] No pepXML found for {d_folder.stem} in {d_parent} or inside {d_folder}")
        return None, None
    if len(pepxmls) > 1:
        print(f"[WARN] Multiple pepXMLs for {d_folder.stem}: {[p.name for p in pepxmls]} — using first.")

    pep_src = pepxmls[0]
    pep_dst = sample_dir / pep_src.name
    if pep_src.resolve() != pep_dst.resolve():
        # move (overwrite if exists)
        pep_src.replace(pep_dst)
    pep_name = pep_dst.name

    # --- 2) find & move matching *_uncalibrated.mzML ---
    mzml_src = pep_src.with_name(f"{pep_src.stem}_uncalibrated.mzML")
    mzml_dst = None
    if not mzml_src.exists():
        # fallback: sometimes Fragger writes just <stem>.mzML
        alt = pep_src.with_name(f"{pep_src.stem}.mzML")
        if alt.exists():
            mzml_src = alt
    if mzml_src.exists():
        mzml_dst = sample_dir / mzml_src.name
        if mzml_src.resolve() != mzml_dst.resolve():
            mzml_src.replace(mzml_dst)
    else:
        print(f"[WARN] No spectra mzML found next to pepXML. Expected {pep_src.stem}_uncalibrated.mzML")

    # --- 3) create alias <pep.stem>.mzML -> *_uncalibrated.mzML ---
    mzml_alias = None
    if mzml_dst and mzml_dst.exists():
        mzml_alias = sample_dir / f"{pep_src.stem}.mzML"
        if not mzml_alias.exists():
            try:
                # relative symlink so the folder stays portable
                os.symlink(mzml_dst.name, mzml_alias)
                print(f"[LINK] {mzml_alias.name} -> {mzml_dst.name}")
            except Exception as e:
                # last resort: duplicate the file (uses disk space)
                shutil.copy2(str(mzml_dst), str(mzml_alias))
                print(f"[COPY] {mzml_alias.name} (symlink failed: {e})")

    # --- 4) handle side files (.mzBIN, .tsv) from the .d parent ---
    for ext in (".mzBIN", ".tsv"):
        src = d_parent / f"{d_folder.stem}{ext}"
        if not src.exists():
            continue
        if side_files == "move":
            dst = sample_dir / src.name
            src.replace(dst)
            print(f"[MOVE] {src.name} -> {dst}")
        elif side_files == "delete":
            src.unlink()
            print(f"[CLEANUP] Removed {src.name}")
        else:
            # skip: leave it where it is
            pass

    print("[DEBUG] Sample dir now contains:", sorted(p.name for p in sample_dir.glob("*")))
    return pep_name, (mzml_alias.name if mzml_alias else None)


def run_msfragger(sample_dir: Path, d_folder: Path, fasta: Path, threads: int):
    """
    sample_dir: working directory for this sample (created)
    d_folder: path to Bruker .d folder
    fasta: FASTA to use
    """
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Creation of fragger.params file
    params_out = sample_dir / "fragger.params"
    write_params_with_fasta(BASE_FRAGGER_PARAMS, params_out, fasta.resolve(), threads)

    # MSFragger reads the .d folder
    cmd = [
        "java", f"-Xmx{JAVA_MEM_GB}G", "-jar",
        str(MSFRAGGER_JAR.resolve()),
        str(params_out),
        str(d_folder.resolve())
    ]
    run(cmd, cwd=sample_dir)

def run_philosopher_chain(sample_dir: Path, fasta: Path, pepxml_names: list[str]):
    """
    pepxml_names: list of filenames (no paths), present in sample_dir.
    """
    if not pepxml_names:
        raise RuntimeError("No pepXML names provided to Philosopher.")

    # init workspace
    run([str(PHILOSOPHER), "workspace", "--init"], cwd=sample_dir)

    # decoy usage
    run([
        str(PHILOSOPHER), "database",
        "--prefix", DECOY_PREFIX,
        "--id", "decoys",
        "--annotate", str(fasta.resolve())
    ], cwd=sample_dir)

    # peptideprophet 
    cmd_pp = [
        str(PHILOSOPHER), "peptideprophet",
        "--expectscore", "--ppm", "--accmass",
        "--decoy", DECOY_PREFIX,
        "--database", str(fasta.resolve()),
    ] + pepxml_names
    run(cmd_pp, cwd=sample_dir)

    # proteinprophet
    interact_pepxml = sorted(p.name for p in sample_dir.glob("interact-*.pep.xml"))
    if not interact_pepxml:
        raise RuntimeError("No interact-*.pep.xml files produced by PeptideProphet.")
    run([str(PHILOSOPHER), "proteinprophet"] + interact_pepxml, cwd=sample_dir)

    prot_xml = sorted(p.name for p in sample_dir.glob("*.prot.xml"))
    if not interact_pepxml:
        raise RuntimeError("No interact-*.pep.xml files found for filter step.")
    if not prot_xml:
        raise RuntimeError("No *.prot.xml file found for filter step.")

    filter_cmd = [
        str(PHILOSOPHER), "filter",
        "--sequential",
        "--tag", DECOY_PREFIX,
        "--razor",
        "--pepxml", *interact_pepxml,   # <---- pepXML
        "--protxml", *prot_xml,         # <---- protXML
        "--psm", "0.01",
        "--pep", "0.01",
        "--prot", "0.01",
    ]

    run(filter_cmd, cwd=sample_dir)
    run([str(PHILOSOPHER), "report"], cwd=sample_dir)

def run_ionquant(spec_dir: str, threads: str, sample_dir: str):
    """
    Runs IonQuant JAR
    """
    psm_path = Path(sample_dir) / "psm.tsv"

    cmd_ionquant = [
        "java", f"-Xmx{JAVA_MEM_GB}G", "-jar", str(IONQUANT_JAR),
        "--specdir", spec_dir,
        "--threads", threads,
        "--ionmobility", "1",
        "--normalization", "0",
        "--mbr", "0",
        "--maxlfq", "0",
        "--psm", str(psm_path)
    ]

    run(cmd_ionquant)

def ms_call(mass_spec_folder: str, fasta: str, result_root: str, threads: int, path_to_d):
    """
    Pipeline to call MS tools
      1. Create work dir
      2. Run MSFragger on .d
      3. Collect outputs into sample_dir
      4. Run Philosopher chain
    """
    d_folder = Path(mass_spec_folder)
    if not d_folder.is_dir() or d_folder.suffix.lower() != ".d":
        return mass_spec_folder, "folder not found or not .d"

    sample_name = d_folder.stem
    sample_dir = Path(result_root) / sample_name
    sample_dir.mkdir(parents=True, exist_ok=True)
    dst = os.path.join(sample_dir, os.path.basename(d_folder))

    #shutil.copytree(d_folder, dst, dirs_exist_ok=True)

    fasta_path = Path(fasta)
    if not fasta_path.exists():
        return mass_spec_folder, f"FASTA not found: {fasta_path}"

    try:
        print(f"[START] {d_folder}")

        td_fasta = sample_dir / f"{Path(fasta).stem}.td.fasta"
        make_target_decoy_fasta(fasta_path, td_fasta, decoy_prefix=DECOY_PREFIX)

        # Run MSFragger with target+decoy FASTA
        run_msfragger(sample_dir, d_folder, td_fasta, threads)
    except Exception as e:
        return mass_spec_folder, f"MSFragger failed: {e}"

    # >>> INSERTED: collect exact pepXML & matching mzML into sample_dir
    pep_name, mzml_alias = collect_msfragger_outputs(sample_dir, d_folder, side_files="move")
    if not pep_name:
        return mass_spec_folder, f"No pepXML found around {d_folder}"
    try:
        # Run Philosopher
        run_philosopher_chain(sample_dir, td_fasta, [pep_name])   # <<< use td_fasta here
    except Exception as e:
        return mass_spec_folder, f"Philosopher failed: {e}"

    try:
        # Run IonQuant
        run_ionquant(spec_dir=str(path_to_d), threads=str(threads), sample_dir=str(sample_dir))
        
    except Exception as e:
        return mass_spec_folder, f"IonQuant failed: {e}"

    print(f"[DONE]  {d_folder}")
    return mass_spec_folder, None

def main():
    pwd = Path(os.getcwd())

    # Input files and directories
    path_to_descriptor = pwd / "scientific_data_processing_files" / "proteome_descript_199_updated.xlsx" # connects the samples and fasta
    foldername = input("Please enter the folder where the .d files are located: ").strip()
    path_to_d = (pwd / foldername).resolve()
    path_to_proteomes = (pwd / "scientific_data_processing_files" / "proteome_fasta_files").resolve()
    zip_locations = (Path(str(path_to_d).replace("PASEF", "_zip_results"))).resolve()
    zip_locations.mkdir(parents=True, exist_ok=True)

    # Parallelization
    total_threads = 128
    threads_per_job = 8
    max_workers = max(1, total_threads // threads_per_job)

    # Loading mapping sheet
    df = pd.read_excel(path_to_descriptor)
    organism_sample_id = df["Filename"].astype(str).tolist()
    proteome_ids = df["Proteome_ID"].astype(str).tolist()

    filenames_fasta_list = os.listdir(path_to_proteomes)
    foldernames_d = [x for x in os.listdir(path_to_d) if x.endswith(".d")]

    # Map proteome ID -> fasta file 
    mapping_id_fasta = {
        id_: next((str(path_to_proteomes / fn) for fn in filenames_fasta_list if id_ in fn), None)
        for id_ in set(proteome_ids)
    }
    # Sample filename -> proteome id
    mapping_sample_proteome_id = dict(zip(organism_sample_id, proteome_ids))

    # Build mapping (.d folder → fasta)
    masspec_to_fasta_dict = {}
    result_root = f"{str(path_to_d)}_msfragger_philosopher_ionquant"
    os.makedirs(result_root, exist_ok=True)

    for dname in foldernames_d:
        curr = dname
        if "_DDA100" in curr:
            curr = curr.split("_DDA100")[0]

        proteome_id = mapping_sample_proteome_id.get(curr)
        fasta_path = mapping_id_fasta.get(proteome_id) if proteome_id else None
        if not fasta_path:
            print(f"[WARN] No FASTA mapped for sample key '{curr}' (from '{dname}'). Skipping.")
            continue

        dpath = str((path_to_d / dname).as_posix())
        masspec_to_fasta_dict[dpath] = fasta_path

    #samples = list(masspec_to_fasta_dict.keys())[:1] # test variant for one sample
    samples = masspec_to_fasta_dict.keys() #NORMAL VERSION
    #samples = [] # rerun failed samples
    print(f"[INFO] Queued samples: {len(samples)}")

    # Run in parallel
    errors = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(ms_call, sample, masspec_to_fasta_dict[sample], result_root, threads_per_job, path_to_d): sample
            for sample in samples
        }
        for fut in as_completed(futures):
            sample, error = fut.result()
            if error:
                print(f"[ERROR] {sample}: {error}")
                errors.append((sample, error))
            else:
                print(f"[OK]    {sample}")

    if errors:
        print("\nSummary of errors:")
        for s, e in errors:
            print(f" - {s}: {e}")
        sys.exit(1)

    print("\nAll MSFragger + Philosopher + IonQuant jobs completed successfully.")


    for result_folder in os.listdir(result_root):
        zip_path = (zip_locations / result_folder).with_suffix(".zip")
        zip_path.parent.mkdir(parents=True, exist_ok=True)
        
        loc_result_folder = f"{result_root}/{result_folder}"
        wanted = {"ion.tsv", "peptide.tsv", "protein.tsv", "psm.tsv"}
        
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for file in os.listdir(loc_result_folder):
                if file in wanted:
                    zf.write(f"{loc_result_folder}/{file}", arcname=file)
        
        print(f"[ZIP] Created {zip_path}")


if __name__ == "__main__":
    main()
