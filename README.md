# daily_tools
Daily tools

<<<<<<< Updated upstream
A collection of daily use tools and scripts to enhance productivity and streamline tasks.

## Features
- Task Automation
- File Management
- Data Processing
- System Monitoring
- Customizable Scripts

## Installation
To install the tools, clone the repository and run the setup script:
```bash
git clone https://github.com/Y-Chao/daily_tools.git
cd daily_tools
python setup.py install
```
=======
## Table of Contents
- **bash** (Bash scripts and tools)
- **cat** (Catalytic related scripts and tools)
- **python** (Python scripts and tools)
- **sci_plots** (Scientific plotting tools)
- **utils** (General utility scripts)

### Connect to the computing node via vscode
Refer to https://github.com/Orion-Zheng/My-HPC-Tunnels/blob/main/launch_tunnel.py

1. Download the "vscode cli" in your remote server, refer [this website](https://code.visualstudio.com/docs/remote/tunnels)
```Bash
cd ~/opt/
curl -Lk 'https://code.visualstudio.com/sha/download?build=stable&os=cli-alpine-x64' --output vscode_cli.tar.gz
tar -xf vscode_cli.tar.gz
rm vscode_cli.tar.gz
```
2. After that, you got a `code` executable file in your `~/opt/` folder.

3. Copy the submit script to the scratch folder, edit the corresponding information (e.g. your conda environment, project name, launch_tunnel.py path, etc.)
```Bash
qsub cpu_tunnel.pbs
```
4. After the job is running, you can use the output file `${RUNTIME_NAME}_tunnel.out`, to follow the instruction to build the tunnel. (e.g. grant access to the server, etc.)

5. Connect to the remote server in the remote vscode application.

### 
>>>>>>> Stashed changes
