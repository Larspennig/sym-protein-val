
# sym-protein-val

## Installation

1. Install Boltz:
   ```bash
   cd modules/boltz
   pip install -e .
   ```

2. Install FoldSeek:
   ```bash
   conda install -c conda-forge -c bioconda foldseek
   ```

3. Install remaining dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. (Optional) Set specific PyTorch CUDA version to work across gtx1080 and a100:
   ```bash
   pip install torch==2.3.1+cu118 -f https://download.pytorch.org/whl/torch_stable.html
   ```


