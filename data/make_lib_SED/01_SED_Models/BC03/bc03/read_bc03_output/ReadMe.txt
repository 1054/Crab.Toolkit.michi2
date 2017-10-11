https://github.com/moustakas/impro/blob/master/pro/synthesis/im_read_bc03.pro
https://raw.githubusercontent.com/moustakas/impro/master/pro/synthesis/im_read_bc03.pro

wget http://www.sdss3.org/svn/repo/idlutils/tags/v5_5_5/pro/misc/struct_trimtags.pro
wget http://www.sdss3.org/svn/repo/idlutils/tags/v5_5_5/pro/misc/struct_addtags.pro
wget http://www.sdss3.org/svn/repo/idlutils/tags/v5_5_8/pro/misc/splog.pro
wget https://idlastro.gsfc.nasa.gov/ftp/pro/structure/mrd_struct.pro


git clone https://github.com/moustakas/impro.git

export bc03_dir="/Users/dzliu/Softwares/BC03/bc03"
idl
bc03 = im_read_bc03(isedfile='../out/out_constant_SFH_2Myr.ised', isedpath='')
