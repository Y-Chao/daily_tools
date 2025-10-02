#!/usr/bin/env python
# -*- encoding: utf-8 -*-

#from __future__ import annotations

__author__ = 'Chao Yang'
__version__=	'1.0'

import os
import sys
import time
import subprocess
from itertools import groupby

ASP2A_Q = "asp2a-login-nus01"

class PBS:
    """
    The job mananeger for PBS pro
    """

    def __init__(self):
        ...

    def monitor(self):
        ...

    def submit(self):
        ...

    def kill(self):
        ...

def split_tasks(tasklist: list, n: int) -> list:
    """
    Split a list of tasks into sublists of size n.
    If the last sublist is smaller than n, it will still be included.
    """
    if not tasklist or n <= 0:
        return []
    split_task = [tasklist[i:i + n] for i in range(0, len(tasklist), n)]
    return split_task

def split_joblist(joblist: str) -> list:
    """
    Split a job list string into a list of individual job names.
    Handles both comma-separated and space-separated lists.
    """
    if not joblist:
        return []
    
    # Replace commas with spaces and split by whitespace
    return [list(g) for k, g in groupby(joblist, key=lambda x: len(x) > 0) if k]
 

def parse_job_info(job_info: str) -> dict:
    """
    Parse the job information string returned by PBS commands.
    Returns a dictionary with job details.
    """
    job_dict = {
        "name": None,
        "id": None,
        "status": None,
        "cores": None,
        "node": None,
        "workdir": None,
        "start": None
    }

    if len(job_info) > 0:
        for line in job_info:
            if "Job Id" in line:
                job_dict["id"] = line.split(":")[1].strip().replace(ASP2A_Q, "")
            elif "Job_Name" in line:
                job_dict["name"] = line.split("=")[1].strip()
            elif "job_state" in line:
                job_dict["status"] = line.split("=")[1].strip()
            elif "Resource_List.ncpus" in line:
                job_dict["cores"] = int(line.split("=")[1].strip())
            elif "Resource_List.nodect" in line:
                job_dict["node"] = int(line.split("=")[1].strip())
            elif "stime" in line:
                t_info = line.split("=")[1].strip()
                t = time.strptime(t_info, "%a %b %d %H:%M:%S %Y")
                t_norm = time.strftime("%Y-%m-%d %H:%M:%S", t)
                job_dict["start"] = t_norm
            elif "exec_host" in line:
                job_dict["node"] = line.split("=")[1].strip()
            elif "PBS_O_WORKDIR" in line:
                line = line.split(",")
                tmp_dict = {}
                for item in line:
                    k = item.split("=")[0]
                    v = item.split("=")[1]
                    tmp_dict.update({k: v})
                workdir = tmp_dict.get("PBS_O_WORKDIR", None)
                if workdir is not None:
                    job_dict["workdir"] = shorten_path(workdir)
            else:
                continue
    return job_dict

def check_jobs(jobname:str=None, jobid:int=None, user:str=None):
    """
    Monitor the status of active jobs in the PBS queue.
    """
    if user is None:
        user = os.getenv("USER")
    try:
        output = subprocess.check_output(["qstat", "-u", user], text=True)
        lines = output.strip().split("\n")[2:]  # Skip header
        active_jobs = [line for line in lines if jobname in line]
        return len(active_jobs)
    except subprocess.CalledProcessError as e:
        print("Error checking jobs:", e)
        return 0
    
def list_all_jobs():
    """
    List all active jobs in the PBS queue.
    """
    cmd = ["qstat", "-f"]
    cmd_a = ["qstat", "-f", "-t"]
    all_jobs = []
    try:
        output = subprocess.check_output(cmd, text=True)
        job_info = output.replace("\n\t", "").split("\n")
        job_list = split_joblist(job_info)
        for job in job_list:
            job_dict = parse_job_info(job)
            if "[]" in job_dict["id"]:
                output_t = subprocess.check_output(cmd_a + [job_dict["id"]], text=True)
                output_t = output_t.replace("\n\t", "").split("\n")
                job_info_t = split_joblist(output_t)
                for job_t in job_info_t:
                    job_dict = parse_job_info(job_t)
                    all_jobs.append(job_dict)
            else:
                all_jobs.append(job_dict)
        return all_jobs
    except subprocess.CalledProcessError as e:
        print("Error listing jobs:", e)
        return []
    
def output_jobs():
    out = list_all_jobs()
    if len(out) == 0:
        print("No jobs")
        return
    headers = list(out[0].keys())
    values = []
    for i, job in enumerate(out):
        values.append(list(job.values()))
    format_print(values, headers)
    return

def monitor_jobs(jobname: str, interval: int = 300):
    """
    Monitor the status of a specific job in the PBS queue.
    """
    if not jobname:
        print("Job name is required for monitoring.")
        return
    from collections import Counter
    run_next = False
    while not run_next:
        counter = Counter()
        all_jobs = list_all_jobs()
        for job in all_jobs:
            if job["name"] == jobname:
                counter[job["status"]] += 1

        if counter["R"] == 0 and counter["Q"] == 0:
            run_next = True
        else:
            print(f"Current job status: {counter}")
            time.sleep(interval)
    print("All jobs completed. Proceeding to the next step.")
    
def submit_job_template(template: str, restart: int = 3):
    """
    Submit a job using a template file.
    """
    for attempt in range(restart):
        try:
            result = subprocess.run(["qsub", template],
                                    check=True,
                                    capture_output=True,
                                    text=True)
            print(f"Job submitted using template: {template}", result.stdout)
            break  # Exit after successful submission
        except subprocess.CalledProcessError as e:
            print("Error submitting job:", e)
            print("Return code:", e.returncode)
            print("Output:", e.output)
            print("Error message:", e.stderr)
            if attempt < restart:
                print(f"Retrying job submission after 10 minutes...")
                time.sleep(600)  # Wait for 10 minutes before retrying
            else:
                raise RuntimeError(f"Failed to submit job after {restart} attempts.")

def shorten_path(path: str):
    home = os.path.expanduser("~")
    short_path = path.replace(home, "~")
    return short_path

def format_print(jobinfo: list, headers: list):
    col_widths = [max(len(str(row[i])) for row in jobinfo) for i in range(len(jobinfo[0]))]
    print("\033[44m" + "  ".join(str(val).center(col_widths[i]) for i, val in enumerate(headers)) + "\033[0m")
    for row in jobinfo:
        print("    ".join(str(val).center(col_widths[i]) for i, val in enumerate(row)))

def batch_sub(num=1):
    """
    Main function to execute the dispatch module.
    This function can be extended to include command-line arguments or other functionalities.
    """
    print("Dispatch module is running...")
    
    # NSCC can run 99 tasks at the same time
    all_list = split_tasks(list(range(num, 5000)), 90)
    list_value = [[l[0], l[-1]] for l in all_list]
    for i, l in enumerate(list_value):
        print(f"Task {i+1}: {l[0]} - {l[1]}")
        if i == 0:
            monitor_jobs("CuSCOH", interval=600)
            submit_job_template("sub_vasp.pbs")
            continue
        os.system(f"sed -i \"8s/{list_value[i-1][0]}/{l[0]}/g\" sub_vasp.pbs")
        os.system(f"sed -i \"8s/{list_value[i-1][1]}/{l[1]}/g\" sub_vasp.pbs")

        # Monitor the job status
        monitor_jobs("CuSCOH", interval=600)
        submit_job_template("sub_vasp.pbs")
    print("All tasks finished.")

def main():
    if len(sys.argv) == 1:
        print("Dispatch module is running...")
        output_jobs()
    else:
        if len(sys.argv) != 3:
            print("Usage: dpdispatch.py [q] [num]")
            sys.exit(1)
        if sys.argv[1] == "q":
            num = int(sys.argv[2])
            batch_sub(num)

if __name__ == '__main__':
   main() 
