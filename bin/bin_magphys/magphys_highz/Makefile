RM = /bin/rm -f
#UTILIB = ./public_utilities.a
#FC =  g77 -ggdb -O0 -ffixed-line-length-132 -fbounds-check
#<DZLIU>#FC = gfortran -std=legacy -ffixed-line-length-none -fbounds-check # need to add "-mcmodel=medium" if we want nfilt>=200
#<DZLIU>#FC = gfortran -std=legacy -ffixed-line-length-none -fbounds-check -mcmodel=medium
FC = gfortran -std=legacy -ffixed-line-length-none -fbounds-check -fno-align-commons -mcmodel=medium

OBJFILES = make_zgrid get_optic_colors_highz get_infrared_colors_highz fit_sed_highz

.f.o:
	$(FC) -c $<

make_zgrid: make_zgrid.f filter.dec $(UTILIB)
	$(FC) make_zgrid.f $(UTILIB) -o make_zgrid

get_optic_colors_highz: get_optic_colors_highz.f filter.dec $(UTILIB)
	$(FC) get_optic_colors_highz.f $(UTILIB) -o get_optic_colors_highz
	
get_infrared_colors_highz: get_infrared_colors_highz.f filter.dec $(UTILIB)
	$(FC) get_infrared_colors_highz.f $(UTILIB) -o get_infrared_colors_highz

fit_sed_highz: fit_sed_highz.f $(UTILIB)
	$(FC) fit_sed_highz.f $(UTILIB) -o fit_sed_highz
	

clean:
	$(RM) public_utilities.a $(OBJFILES)

all:
#	make public_utilities.a
	make make_zgrid get_optic_colors_highz get_infrared_colors_highz fit_sed_highz

#---------------------------------------------------------------------------
#commands to build library public_utilities.a
#SRCFILES = make_zgrid.f get_optic_colors_highz.f get_infrared_colors_highz.f fit_sed_highz.f

#OBJFILES = make_zgrid.o get_optic_colors_highz.o get_infrared_colors_highz.o fit_sed_highz.o

#public_utilities.a: $(OBJFILES)
#	ar ruv public_utilities.a $?
#	ranlib public_utilities.a
