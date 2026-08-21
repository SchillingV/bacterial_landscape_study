import os
import sys
import shutil
import pandas as pd
import subprocess as sub
from concurrent.futures import ProcessPoolExecutor, as_completed

### DIA-NN Commands ###

def command_1(rel_path, mass_spec_folder, mass_spec_folder_no_space_no_d, fasta, threads):
    command = [
    "diann",
    "--f", mass_spec_folder,
    "--lib", "",
    "--threads", str(threads),
    "--verbose", "1",
    "--out", f"{rel_path}/{mass_spec_folder_no_space_no_d}_report.tsv",
    "--qvalue", "0.01",
    "--out-lib", f"{rel_path}/{mass_spec_folder_no_space_no_d}-lib.tsv",
    "--gen-spec-lib",
    "--predictor",
    "--fasta", fasta,
    "--fasta-search",
    "--min-fr-mz", "100",
    "--max-fr-mz", "1700",
    "--met-excision",
    "--cut", "K*,R*,!*P",
    "--missed-cleavages", "1",
    "--min-pep-len", "7",
    "--max-pep-len", "30",
    "--min-pr-mz", "350",
    "--max-pr-mz", "1150",
    "--min-pr-charge", "2",
    "--max-pr-charge", "4",
    "--unimod4",
    "--individual-mass-acc",
    "--individual-windows",
    "--relaxed-prot-inf",
    "--smart-profiling",
    "--pg-level", "1",
    "--peak-center",
    "--no-ifs-removal"]

    try:
        sub.run(command, check=True)
    except sub.CalledProcessError as e:
        print(f"Error occurred while running the first command: {e}")

def command_2(rel_path, mass_spec_folder_no_space_no_d, threads):
    command2 = [
    "diann",
    "--lib", f"{rel_path}/{mass_spec_folder_no_space_no_d}-lib.predicted.speclib",
    "--threads", str(threads),
    "--verbose", "1",
    "--out", f"{rel_path}/{mass_spec_folder_no_space_no_d}.tsv",
    "--qvalue", "0.01",
    "--matrices",
    "--out-lib", f"{rel_path}/{mass_spec_folder_no_space_no_d}-lib.predicted.tsv",
    "--gen-spec-lib",
    "--relaxed-prot-inf",
    "--smart-profiling",
    "--peak-center",
    "--no-ifs-removal"]

    try:
        sub.run(command2, check=True)
    except sub.CalledProcessError as e:
        print(f"Error occurred while running the second command: {e}")

def diann_call(mass_spec_folder, fasta, result_folder, threads):
    # Check folder exists
    if not os.path.isdir(mass_spec_folder):
        return mass_spec_folder, "folder not found"

    # build “clean” name & create output dir
    mass_spec_folder_no_space_no_d = mass_spec_folder.replace(".d", "").replace("BacLandscape_", "").split("/")[-1]
    rel_path = f"{result_folder}/{mass_spec_folder_no_space_no_d}"
    os.makedirs(f"{result_folder}/{mass_spec_folder_no_space_no_d}", exist_ok=True)

    # run DIA‑NN command 1
    try:
        command_1(rel_path, mass_spec_folder, mass_spec_folder_no_space_no_d, fasta, threads)
    except Exception as e:
        return mass_spec_folder, f"command_1 failed: {e}"

    # run DIA‑NN command 2 to save spec library
    try:
        command_2(rel_path, mass_spec_folder_no_space_no_d, threads)
    except Exception as e:
        return mass_spec_folder, f"command_2 failed: {e}"

    return mass_spec_folder, None

def zip_folders():
    pass

def main():

    # Get locations for the excel sheet and fast files 

    pwd = os.getcwd()
    path_to_descriptor = os.path.join(pwd, "proteome_descript_199.xlsx") # sample to fasta info
    foldername = input("Please enter the folder where the .d files are located: ")
    path_to_d = os.path.join(pwd,foldername)
    path_to_proteomes = os.path.join(pwd,"/proteome_fasta_files/")
    zip_locations = os.path.join(path_to_d.replace("PASEF", "zip_files"))


    # threading strategy and job parallelization - adjust this for reanalysis of the data depending on your system
    total_threads    = 128
    threads_per_job  = 16
    max_workers      = total_threads // threads_per_job  # 24 or 32

    # Get list of mass spec folders and fasta files; file operations can be probably removed or adjusted depending on your paths and system 

    df = pd.read_excel(path_to_descriptor)

    organism_sample_id = df["Filename"].to_list()
    proteome_ids = df["Proteome_ID"].to_list()
    filenames_fasta_list = os.listdir(path_to_proteomes)
    foldernames_d = os.listdir(path_to_d)

    mapping_id_fasta = {id_: filename for id_ in list(set(proteome_ids)) for filename in filenames_fasta_list if id_ in filename}
    mapping_sample_proteome_id = dict(zip(organism_sample_id, proteome_ids))

    masspec_to_fasta_dict = {}
    result_folder = f"{path_to_d}_diann_res"
    os.makedirs(result_folder, exist_ok=True)
    os.makedirs(zip_locations, exist_ok=True)

    for foldername in foldernames_d:
        curr_sample = foldername.split("PASEF_")[1]
        curr_sample = curr_sample.split("-")[0]
        curr_sample = curr_sample.rsplit("_", 1)[0]
        folderpath = os.path.join(path_to_d, foldername).replace("\\", "/")
        masspec_to_fasta_dict[folderpath] = os.path.join(path_to_proteomes.replace("\\", "/"), mapping_id_fasta[mapping_sample_proteome_id[curr_sample]])
    #print(mapping, len(mapping))

    samples = list(masspec_to_fasta_dict.keys())

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # submit one future per sample
        futures = {
            executor.submit(diann_call, sample, masspec_to_fasta_dict[sample], result_folder, threads_per_job): sample
            for sample in samples
        }
        for fut in as_completed(futures):
            sample, error = fut.result()
            if error:
                print(f"[ERROR] {sample}: {error}")
            else:
                print(f"[DONE]  {sample}")

    print("All DIA-NN pipelines completed.")

    ### Copy .quant files to result folders and zip .d files for completion

    for folder in foldernames_d:

        ### copy quant files to correct location and keep naming convention within the folder
        quant_loc_old = os.path.join(pwd, path_to_d, f"{folder}.quant")
        quant_loc_new = os.path.join(pwd, result_folder, f"{folder.replace("BacLandscape_", "")[:-2]}", f"{folder.replace("BacLandscape_", "")[:-2]}.quant")

        shutil.copyfile(quant_loc_old, quant_loc_new)

        ### zip files for later upload to PRIDE ###
        loc_d_folder = os.path.join(pwd, path_to_d, folder)
        zip_base = os.path.join(pwd, zip_locations, f"{folder.replace("BacLandscape_", "")}")
        shutil.make_archive(zip_base, "zip", loc_d_folder)
        #shutil.copyfile(os.path.join(zip_base, "zip"), os.path.join(zip_locations, zip_base, f"{folder.replace("BacLandscape_", "")}", "zip"))
    
    print("ALL IS FINISHED")
    

if __name__ == "__main__":
    main()
