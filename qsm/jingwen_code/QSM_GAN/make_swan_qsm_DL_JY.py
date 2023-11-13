'''
This is the script to calculate QSM from tissue phase using DL-based QSM algorithm.

Example usage:
    > python make_swan_qsm_DL.py /data/7T_hunt/b4603/t12772/swan_qsm/t12772_tissue_phase.idf

Eason Chen
5/23/2019
'''
import glob, os
import json
from predict import *
import sys
import time

start_time = time.time()
args = sys.argv
assert len(args) == 2, "Please specify ONE tissue phase filepath"
idf = args[1]

print('->Checking local field file type...')
if '.nii.gz' in idf:
    print('--->Nifti file already exist')
else:
    print('->Converting idf to nifti...')
    try:
        output = os.path.dirname(idf)
        output_root = os.path.basename(idf).split('.')[0]
        os.system(f'nifti_file_convert --input {idf} --output {output} --output_root {output_root}')
    except:
        raise NameError('Convert idf to nifti failed!')

# select network weights and config file to load
ckpt_path = '../exp/best_model/WGAN_i64o48/net_best.pt'
cfg_path  = '../exp/best_model/WGAN_i64o48/config.json'

print('->Loading model weights...')
with open(os.path.join(cfg_path), 'r') as f:
    cfg = json.load(f)
device = cfg['device']

net = UNet3D(**cfg['netG_args']).to(device)
ckpt = torch.load(ckpt_path)
net.load_state_dict(ckpt['netG'])

print('->Calculating QSM and saving to nifti...')
nifti_files = [idf.split('.')[0] + '.nii.gz']
qp = QsmPredictHD(cfg, nifti_files)
qp.run_net(net)
qp.save_predicts()

end_time = time.time()

print(f'Total time elapsed: {end_time - start_time:.2}s')
