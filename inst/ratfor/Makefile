
RATFOR = /home/adrian/bin/ratfor77
#RATFOR = /usr/local/bin/ratfor
CPP = /usr/bin/cpp 

##########################################################
# Sources actually written by humans:

RAT_SRC = dppll.r inxypOld.r 

C_DOMINIC = dinfty.c dwpure.c 

C_MISC = raster.h areadiff.c closepair.c connect.c corrections.c \
       discarea.c distances.c distmapbin.c distseg.c \
       exactdist.c  exactPdist.c  \
       massdisthack.c poly2im.c trigraf.c utils.c xyseg.c 

C_MH = methas.h dist2.h areaint.c badgey.c dgs.c \
	diggra.c dist2.c fexitc.c getcif.c geyer.c \
	lookup.c methas.c stfcr.c \
	straush.c straushm.c strauss.c straussm.c

C_KEST = Kloop.h Kborder.c

C_SRC = $(C_DOMINIC) $(C_MISC) $(C_MH) $(C_KEST)

CC_SRC = PerfectStrauss.cc

HUMAN = $(RAT_SRC) $(C_SRC) $(CC_SRC)  Makefile

##########################################################
# Source to be generated automatically:

RAT_FOR = dppll.f inxypOld.f 

GENERATED = $(RAT_FOR) 


######################################################
###########  TARGETS ################################

target:	$(GENERATED)
	@echo -- Done -------

tar:
	tar cvf src.tar $(HUMAN)

clean:
	rm $(GENERATED)
	-rm src.tar

#######################################################
#########    RULES   ##################################


.r.f:
	$(RATFOR) -o $@ $?

