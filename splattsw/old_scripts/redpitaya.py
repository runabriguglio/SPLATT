import os

def set_frequency(freq):
    cmd=r'sshpass -p root ssh root@193.206.155.193 /mnt/scripts/obb.sh %f ' % freq
    retval = os.system(cmd)
