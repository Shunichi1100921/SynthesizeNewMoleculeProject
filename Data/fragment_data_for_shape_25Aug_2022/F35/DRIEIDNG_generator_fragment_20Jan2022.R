## R script to generate LAMMPS input file for the DRIEDING force field
## This version generates a box filled with molecules in random orientations.
##

## Update on 29 Oct 2021: Designed for input from Molecule fragmenter code (molecule_fragmenter_code_21Oct2021.R). This version of the code does not use PDB file as input.

## Based on DRIEDING_generator.R (3-7-19, 2-4-20, 15-9-20, 29-10-21, 12-1-22)

## Update on 12 Jan 2022: Fixed error which caused Cl atom to be associated with an angle.

## Update on 20 Jan 2022: Fixed another error associated with Cl, Br, and I (atomic masses were not properly included)

## By Daniel Packwood (start: 23-12-20)

rm(list = ls())

## Geometry and vdW parameters (need to input by hand)

## Basic input / output

at.types.file <- "attypes.dat" # Name of file which contains the atom types
bond.file <- "bond.dat" # Name of the file which contains the matrix of bonds
coord.file <- "coord.dat" # Name of the file which contains the coordinate matrix

xlo <- 0; xhi <- 15.0 # Simulation box dimensions (angstroms)
ylo <- 0; yhi <- 15.0 # Simulation box dimensions (angstroms)
zlo <- 0; zhi <- 15.0 # Simulation box dimensions (angstroms)

xy <- 0.0; # Incase a non-cubic simulation box is wanted
xz <- 0.0; yz <- 0.0;

dat.output <- "frag5.dat" # Name of the .dat file to output
in.output <- "frag5.in" # Name of the .in file to output

autowrite.dat <- TRUE # If TRUE, the .dat file is written automatically upon execution
autowrite.in <- TRUE # If TRUE, the .in file is written automatically upon execution

nMol <- 1 # Number of molecules to put into simulation box (must be fixed at present...)
rMin <- 1 # Minimum distance (angstroms) between any pair of molecules placed in the simulation box (for generating the initial condition)

nWat <- 500 # Number of water molecules to include

############################################################################################
############### PART 1 - Determine the atoms, bonds, torsions, and impropers ###############
############################################################################################

## Read and process atom types

AtTypes <- readLines(at.types.file)

nAt <- length(AtTypes) # Number of atoms in the molecule.

AtTypesU <- unique(AtTypes)
AtIndex <- 1:length(AtTypes)

if(nWat > 0 & !"H__HB"%in%AtTypesU)
	AtTypesU <- c(AtTypesU, "H__HB")
	
if(nWat > 0 & !"O_3"%in%AtTypesU)
	AtTypesU <- c(AtTypesU, "O_3")
	
## Assign atom masses

massT <- c()
for(k in 1:length(AtTypesU))
	{if(grepl("C", AtTypesU[k]))
		{massT <- c(massT, 12.01)

		 if(grepl("Cl", AtTypesU[k]))
			massT[length(massT)] <- 35.45
		}

	 if(grepl("N", AtTypesU[k]))
	 	massT <- c(massT, 14.01)

	 if(grepl("O", AtTypesU[k]))
	 	massT <- c(massT, 16.00)

	 if(grepl("H", AtTypesU[k]))
		massT <- c(massT, 1.01)
		
	 if(grepl("F", AtTypesU[k]))
		massT <- c(massT, 19.00)

	 if(grepl("S", AtTypesU[k]))
		massT <- c(massT, 32.06)

	 if(grepl("Br", AtTypesU[k]))
		massT <- c(massT, 79.90)

	 if(grepl("I", AtTypesU[k]))
		massT <- c(massT, 126.91)

	 if(grepl("P", AtTypesU[k]))
		massT <- c(massT, 30.97)
	}

## Read and process geometric valance parameters

make_gvp <- function() # Generate matrix of geometric valance parameters
	{gvpm <- array(0,c(0,3))
	
	 if("H_"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("H_", "0.330", "180.0"))
		
	 if("H__HB"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("H__HB", "0.330", "180.0"))
	 
	 if("C_3"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("C_3", "0.770", "109.471"))
	 
	 if("C_R"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("C_R", "0.700", "120.0"))

	 if("C_2"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("C_2", "0.670", "120.0"))
		
	 if("C_1"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("C_1", "0.602", "180.0"))
			
	 if("N_3"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("N_3", "0.702", "106.7"))
	
	 if("N_R"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("N_R", "0.650", "120.0"))

	 if("N_2"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("N_2", "0.615", "120.0"))

	 if("N_1"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("N_1", "0.556", "180.0"))

	 if("O_3"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("O_3", "0.660", "104.51"))
		
	 if("O_R"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("O_R", "0.660", "120.0"))
		
	 if("O_2"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("O_2", "0.560", "120.0"))
		
	 if("O_1"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("O_1", "0.528", "180.0"))

	 if("F_"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("F_", "0.611", "180.0"))

	 if("P_3"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("P_3", "0.890", "92.3"))

	 if("S_3"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("S_3", "1.040", "92.1"))
		
	 if("Cl_"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("Cl_", "0.997", "180.0"))

	 if("Br_"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("Br_", "1.167", "180.0"))

	 if("I_"%in%AtTypesU)
		gvpm <- rbind(gvpm, c("I_", "1.360", "180.0"))

	 gvpm	
	}

gvp.dat <- make_gvp()

if(nWat > 0 & !"H__HB"%in%gvp.dat[,1])
	gvp.dat <- rbind(gvp.dat, c("H__HB", "0.330", "180.0"))

if(nWat > 0 & !"O_3"%in%gvp.dat[,1])
	gvp.dat <- rbind(gvp.dat, c("O_3", "0.660", "104.51"))
		
## Read and process the van der Waals parameters

get_vdwp <- function(At) # Get vdW parameters for atom type At
	{at <- strsplit(At,"_")[[1]][1]
		
	 #if(at == "H")
		#vout <- t(c(At, "3.195", "0.0152")) 

	 if(At == "H_")
	       vout <- t(c(At, "3.195", "0.0152"))

	 if(At == "H__HB")
	       vout <- t(c(At, "3.193", "0.0001"))

	 if(at == "C")
		vout <- t(c(At, "3.8983", "0.0951"))

	 if(at == "N")
		vout <- t(c(At, "3.6621", "0.0774"))

	 if(at == "O")
		vout <- t(c(At, "3.4046", "0.0957"))

	 if(at == "F")
		vout <- t(c(At, "3.4720", "0.0725")) 
		
	 if(at == "P")
		vout <- t(c(At, "4.15", "0.32")) 
		
	 if(at == "S")
		vout <- t(c(At, "4.03", "0.344")) 
		
	 if(at == "Cl")
		vout <- t(c(At, "3.95", "0.2833"))
		
	 if(at == "Br")
		vout <- t(c(At, "3.95", "0.37")) 
		
	 if(at == "I")
		vout <- t(c(At, "4.15", "0.51")) 
		
	 vout	
	}

make_vdw <- function() # Generate matrix of van der Waals parameters
	{vdwm <- array(0,c(0,3))
	
	 for(k in AtTypesU)
		vdwm <- rbind(vdwm, get_vdwp(k))
		
	 vdwm	
	}

vdW.dat <- make_vdw()

# Obtain coordinate information

Zcoord <- read.table(coord.file)
Zcoord <- as.matrix(Zcoord)

# Generate initial condition: Create nMol copies of Zcoord (one for each molecule), subject to no-overlap constraints
# Note on 15-9-20: Need to add a routine to rotate the molecules

LatVec <- array(0, c(3,3))
LatVec[1,1] <- xhi - xlo; 
LatVec[2,1] <- xy; LatVec[2,2] <- yhi - ylo; LatVec[2,3] <- 0;
LatVec[3,1] <- xz; LatVec[3,2] <- yz; LatVec[3,3] <- zhi - zlo;

# Rotation matrices

rot_mat <- function(tz1, ty1, tzd1) # Make rotation matrix with Euler angles tz1, ty1, tzd1
	{
	 Rz1 <- array(0, c(3,3)); Ry1 <- array(0, c(3,3)); Rzd1 <- array(0,c(3,3))

	 Rz1[1,1] <- cos(tz1); Rz1[1,2] <- -sin(tz1); Rz1[1,3] <- 0;
	 Rz1[2,1] <- sin(tz1); Rz1[2,2] <- cos(tz1); Rz1[2,3] <- 0;
	 Rz1[3,1] <- 0; Rz1[3,2] <- 0; Rz1[3,3] <- 1;

	 Ry1[1,1] <- cos(ty1); Ry1[1,2] <- 0; Ry1[1,3] <- sin(ty1);
	 Ry1[2,1] <- 0; Ry1[2,2] <- 1; Ry1[2,3] <- 0;
	 Ry1[3,1] <- -sin(ty1); Ry1[3,2] <- 0; Ry1[3,3] <- cos(ty1);

	 Rzd1[1,1] <- cos(tzd1); Rzd1[1,2] <- -sin(tzd1); Rzd1[1,3] <- 0;
	 Rzd1[2,1] <- sin(tzd1); Rzd1[2,2] <- cos(tzd1); Rzd1[2,3] <- 0;
	 Rzd1[3,1] <- 0; Rzd1[3,2] <- 0; Rzd1[3,3] <- 1;

	 Rot1 <- Rz1; Rot1 <- Ry1%*%Rot1; Rot1 <- Rzd1%*%Rot1;
	 
	 Rot1
	}

# Functions for movement

get_com <- function(dx, dy, dz) # Convert direct COM coordinates (dx, dy, dz) into real coordinates
	{
	 Dcom <- as.matrix(c(dx, dy, dz))
	 Rcom  <- LatVec%*%Dcom
	
	 c(Rcom)
	}
		
rotate_move <- function(comx, comy, comz, tz1, ty1, tzd1) # Rotate and move a molecule so that: center-of-mass is (comx, comy, comz), and Euler rotation angles are (tz1, ty1, tzd1)
	{Rotm <- rot_mat(tz1, ty1, tzd1)
	
	COM <- get_com(comx, comy, comz)
				
	Zout <- Zcoord%*%Rotm;
	
	# Move molecule to new COM position
	
	com_n <- apply(Zout, MARGIN = 2, mean);
	
	Zout[,1] <- Zout[,1] - com_n[1] + COM[1];
	Zout[,2] <- Zout[,2] - com_n[2] + COM[2];
	Zout[,3] <- Zout[,3] - com_n[3] + COM[3];	
	
	Zout	
	}	

mat.dist <- function(Z1, Z2) # Find minimum distance between two matrices Z1 and Z2
	{Z12 <- rbind(Z1, Z2)

	 dZ12 <- dist(Z12)
	 dZ12 <- as.matrix(dZ12)

	 dZ12 <- dZ12[1:nrow(Z1), (nrow(Z1) + 1):(nrow(Z1) + nrow(Z2))]

	 min(dZ12)
	}

make.ic <- function() # Place the molecules in the box
	{Zlist <- list()	
	 OK = FALSE
	 
	 while(!OK)
		{OK = TRUE
		
		 for(k in 1:nMol)
			{cx = runif(1,0,1)	
			 cy = runif(1,0,1)
			 cz = runif(1,0,1)
			 
			 tz1 <- runif(1,0,2*pi)
			 ty1 <- runif(1,0,2*pi)
			 tzd1 <- runif(1,0,2*pi)
			
			 Zlist[[k]] <- rotate_move(cx, cy, cz, tz1, ty1, tzd1)			
			}		
			
		 if(nMol > 1)
		{
		 for(k in 1:(nMol - 1))
			{
			for(j in (k + 1):nMol)
				{dkj <- mat.dist(Zlist[[k]], Zlist[[j]])
				
				 if(dkj < rMin)
					{OK = FALSE
					 break		
					}
				}
				
			if(!OK)
				break
			}
		}
		
		}
		
	Zlist
	}

ZcoordL <- make.ic() #make.ic()
print("Molecule coordinates made")

# Incorporate water molecules

# Create water coordinates

get.h2o <- function()
	{Zout <- array(0, c(3,3))

	 Zout[1,1] <- -0.96; Zout[1,2] <- 0; Zout[1,3] <- 0;
	 Zout[2,1] <- 0; Zout[2,2] <- 0; Zout[2,3] <- 0;
	 Zout[3,1] <- 0.24; Zout[3,2] <- -0.93; Zout[3,3] <- 0;
	 
	 Zout
	}	

Wcoord <- get.h2o()

move.com <- function(Zcoord, Rcom) # Move Zcoord so that molecule center-of-mass is at point Rcom
	{com1 <- c(mean(Zcoord[,1]), mean(Zcoord[,2]), mean(Zcoord[,3]))

	 Zcoord2 <- Zcoord;
	 Zcoord2[,1] <- Zcoord[,1] + Rcom[1] - com1[1]
	 Zcoord2[,2] <- Zcoord[,2] + Rcom[2] - com1[2]
	 Zcoord2[,3] <- Zcoord[,3] + Rcom[3] - com1[3]

	 Zcoord2
	}

make.ic.water <- function() # Make initial conditions for water molecules. Overlaps are not accounted for here.
	{WcoordL <- list()

	 for(k in 1:nWat)
		{Wtemp <- move.com(Wcoord, c(runif(1,xlo,xhi), runif(1,ylo,yhi), runif(1,zlo,zhi)))
		
		 WcoordL[[k]] <- Wtemp
		}

	 WcoordL	 	 
	}
	
if(nWat > 0)   
	{WcoordL <- make.ic.water()
	 print("Initial water coordinates made")
	}	

# Obtain bonding information

bondm <- read.table(bond.file)
bond.mat <- as.matrix(bondm)
	
get.bond.mat.h2o <- function() # Generate the bonding matrix for water
	{wbmat <- array(0,c(2,2))

	 wbmat[1,1] <- 1; wbmat[1,2] <- 2;
	 wbmat[2,1] <- 2; wbmat[2,2] <- 3;

	 wbmat
	}

if(nWat > 0)
	wbond.mat <- get.bond.mat.h2o()

nBond <- nrow(bond.mat)

if(nWat > 0)
	nBondW <- nrow(wbond.mat)

## Determine bond types and parameters

bond.dat <- array("", c(nBond*nMol, 6)) # Data for each bond. first and second columns are for the two atoms in the bond, third column is bond type (single, double, or resonance), fourth column is equilibrium bond length, third column is bond type, and the fifth column is force constant
bond.dat.single <- array("", c(nBond, 6)) # As for bond.dat, but for a single molecule

get.bl <- function(t1,t2) # Get the bond length for a bond between atoms of the type t1 and t2
	{wh1 <- which(gvp.dat[,1] == t1)
	 wh2 <- which(gvp.dat[,1] == t2)

	 bl12 <- as.numeric(gvp.dat[wh1,2]) + as.numeric(gvp.dat[wh2,2]) - 0.01 # 0.01 is the DRIEDING 'delta' parameter

	 bl12
	}

get.sdr <- function(t1,t2) # Determine if bond is single, double, or resonance
	{# Get last bond indices

	 t1_last <- strsplit(t1, "_")[[1]]
	 t1_last <- t1_last[length(t1_last)]

	 t2_last <- strsplit(t2, "_")[[1]]
	 t2_last <- t2_last[length(t2_last)]

	 btype <- "single"

	 # Condition for double (not connected to an sp3 atom)

	 if(!"3"%in%c(t1_last, t2_last))
		btype <- "double"

	 # Condition for resonance (both atoms are type "R")

	 if(t1_last == "R" & t2_last == "R")
		btype <- "resonance"

	 # Condition for single bond involving H

	 if(c("H")%in%c(t1_last, t2_last) | c("HB")%in%c(t1_last, t2_last))
		btype <- "single"

	 # Condition for single bond between sp2 C and resonance C

	 if(t1_last == "2" & t2_last == "R")
		btype <- "single"

	 if(t2_last == "2" & t1_last == "R")
		btype <- "single"

	 # Condition for C(sp2)-N(sp2) single bond. Will need to re-consider this if C=N double bonds are ever present in the compound

	 #if(t1 == "C_2" & t2 == "N_2")
		#btype <- "single"

	 #if(t1 == "N_2" & t2 == "C_2")
		#btype <- "single"	 

	 # Condition for N(sp2)-N(sp2) single bond. Will need to re-consider this if N=N double bonds are ever present in the compound

	 if(t1 == "N_2" & t2 == "N_2")
		btype <- "single"


	 btype
	}

count <- 1

for(k in 1:nrow(bond.mat))
	{# Determine bond type

	 i1 <- bond.mat[k,1]
	 i2 <- bond.mat[k,2]

	 t1 <- AtTypes[i1]
	 t2 <- AtTypes[i2]

	 wh1 <- which(AtTypesU == t1)
	 wh2 <- which(AtTypesU == t2)

	 o1 <- which.min(c(wh1,wh2))
	 o2 <- which.max(c(wh1,wh2))

	 t12 <- paste(c(t1,t2)[c(o1,o2)], collapse = "*")

	 # Determine equilibrium bond length

	 rIJ = get.bl(t1,t2)

	 # Determine whether the bond is single or double

	 btype <- get.sdr(t1,t2)

	 # Force constant

	 KIJ <- 700 # Default is 700 (kcal/mol)/Angstrom
		
	 if(btype == "double") 
		KIJ = 2*700

	 if(btype == "resonance")
		KIJ = 1.5*700 # Is this the correct bond order?

	 # Input values to matrix, one for each molecule

	 bond.dat.single[k,] <- c(i1,i2, t12, btype, rIJ, KIJ)

	 for(j in 1:nMol)
		{bond.dat[count,] <- c(i1 + (j-1)*nrow(Zcoord), i2 + (j-1)*nrow(Zcoord), t12, btype, rIJ, KIJ)
		 count <- count + 1
		}
	}
	
if(nWat > 0)
{
nAt_start <- nAt*nMol + 1 # Atom numbering to start counting for H20 molecules
count <- nAt_start

bond.dat.h2o <- array(0,c(0,6))
AddL <- array(0,c(2,6))

for(k in 1:nWat)
	{
	 AddL[1,1] <- count; AddL[1,2] <- count + 1; AddL[1,3] <- "O_3*H__HB"; AddL[1,4] <- "single"; AddL[1,5] <-  get.bl("O_3","H__HB"); AddL[1,6] <- 700
	 AddL[2,1] <- count + 1; AddL[2,2] <- count + 2; AddL[2,3] <- "O_3*H__HB"; AddL[2,4] <- "single"; AddL[2,5] <-  get.bl("O_3","H__HB"); AddL[2,6] <- 700

	 count <- count + 3

	 bond.dat.h2o <- rbind(bond.dat.h2o, AddL)
	}

#Don't need thiR

#bond.dat <- rbind(bond.dat, bond.dat.h2o)
}

## Determine angle types and parameters
## We use harmonic angles in this set-up (eq. 11 of Mayo et al), as this is easier to implement in LAMMPS

column.flip <- function(mat) # Switch the places of the first two columns of matrix mat
	{mato <- mat
	 mato[,2] <- mat[,1]
	 mato[,1] <- mat[,2]

	 mato
	}

get.partners <- function(k) # Determine the angles that atom k is involved in
	{wh1 <- which(bond.mat[,1] == k)
	 wh2 <- which(bond.mat[,2] == k)

	 bondk1 <- bond.mat[wh1,]
	 bondk2 <- bond.mat[wh2,]

	 if(length(bondk1) == 2)
		bondk1 <- t(bondk1)

	 if(length(bondk2) == 2)
		bondk2 <- t(bondk2)

	 bondp <- c(bondk1[,2], bondk2[,1]) # Atoms bonded to k

	 bondp
	}

get.angles <- function(k) # Get the angle data for angles centered on k
	{bondp <- get.partners(k)

	 bondp_combn <- t(combn(bondp,2))
	 angle_k <- cbind(bondp_combn[,1],k,bondp_combn[,2])

	 # Set angle type

	 atype <- AtTypes[AtIndex[k]]

	 # Set equilibrium angle

	 wht <- which(gvp.dat[,1] == atype)
	 eqang <- gvp.dat[wht,3]

	 # Set force constant

	 Kang <- 100 # (kcal/mol)/rad^2

	 # Output

	 cbind(atype,Kang,eqang,angle_k)
	}

exclude_cases <- c("H_", "H__HB", "O_2", "F_", "Cl_", "Br_", "I_") # Atoms which cannot contribute a bond angle (since they are bonded to one partner only)

ang.dat <- array(0,c(0,6))

for(k in 1:nAt) # Compile the data on bond angles
	{ktype <- AtTypes[AtIndex[k]]

	 if(!ktype%in%exclude_cases)
	 	{#ang.dat <- rbind(ang.dat, get.angles(k))

		 ang.dat.add <- get.angles(k)
		 
		 for(j in 1:nMol)
		 	{ang.dat.k <- ang.dat.add

			 ang.dat.k[,4] <- as.numeric(ang.dat.add[,4]) + (j-1)*nAt
			 ang.dat.k[,5] <- as.numeric(ang.dat.add[,5]) + (j-1)*nAt
			 ang.dat.k[,6] <- as.numeric(ang.dat.add[,6]) + (j-1)*nAt

			 ang.dat <- rbind(ang.dat, ang.dat.k)
		 	}
		}
	}

if(nWat > 0)
{
ang.dat.h2o <- array(0,c(nWat,6)) # Angle data for H2O

count <- 1

for(k in 1:nWat)
	{at1 <- bond.dat.h2o[count,1];
 	 at2 <- bond.dat.h2o[count,2];
	 at3 <- bond.dat.h2o[(count+1),2];

	 ang.dat.h2o[k,1] <- "O_3";  ang.dat.h2o[k,2] <- "100"; ang.dat.h2o[k,3] <- "104.51"; 
	 ang.dat.h2o[k,4] <- at1; ang.dat.h2o[k,5] <- at2; ang.dat.h2o[k,6] <- at3;

	 count <- count + 2
	}

# Don't need this either	
# ang.dat <- rbind(ang.dat, ang.dat.h2o)
}

## Determine dihedral angle data

dihedral_exists <- function(h) # Check whether dihedral angle exists for a bond
	{htype <- bond.dat.single[h,3]

	 dexists <- TRUE

	 if(grepl("H", htype))
		dexists <- FALSE

	 if(grepl("O_2", htype))
		dexists <- FALSE
		
	 if(grepl("F", htype))
		dexists <- FALSE

	 dexists
	}

wh.dihedral <- sapply(1:nrow(bond.dat.single), dihedral_exists, USE.NAMES = FALSE)
wh.dihedral <- c(1:nrow(bond.dat.single))[wh.dihedral]

dindex.get <- function(h) # Get the bonds connected to bond h (h is a row of bond.dat.single)
	{bondh <- bond.mat[h,]

	 i1 <- bondh[1]
	 i2 <- bondh[2]

	 # Find all bonds involving atom i1

	 wh1 <- which(bond.mat[,1] == i1)
	 wh2 <- which(bond.mat[,2] == i1)

	 wh12 <- c(wh1,wh2)
	 wh12 <- wh12[!wh12==h]

	 # Put these bonds into format that is easy for creating dihedral data

	 bond_12 <- bond.mat[wh12,]

	 if(length(bond_12) == 2)
		bond_12 <- t(bond_12)	 

	 for(i in 1:nrow(bond_12))
		if(bond_12[i,1] == i1)
			bond_12[i,] <- bond_12[i,c(2,1)]

	 # Find all bonds involving atom i2

	 whA <- which(bond.mat[,1] == i2)
	 whB <- which(bond.mat[,2] == i2)

	 whAB <- c(whA,whB)
	 whAB <- whAB[!whAB==h]

	 # Put these bonds into format that is easy for creating dihedral data

	 bond_AB <- bond.mat[whAB,]

	 if(length(bond_AB) == 2)
		bond_AB <- t(bond_AB)	 

	 for(i in 1:nrow(bond_AB))
		if(bond_AB[i,2] == i2)
			bond_AB[i,] <- bond_AB[i,c(2,1)]

	 # Put into dihedral format

	 b1 <- 1:nrow(bond_12)
	 b2 <- 1:nrow(bond_AB)

	 b12 <- expand.grid(b1,b2)

	 di.dat <- array(0,c(0,4))

	 for(i in 1:nrow(b12))
		{bond1 <- bond_12[b12[i,1],]
	 	 bondA <- bond_AB[b12[i,2],]

		 addi <- t(c(bond1,bondA))

		 di.dat <- rbind(di.dat, addi)
		}

	 # Output

	 di.dat
	}


get.dihedral.order <- function(h) # Get the order of the bond in the middle of the torsion
	{a1 <- dihedral_array[h,2]; a2 <- dihedral_array[h,3]

	 wha1 <- which(bond.dat.single[,1] == as.character(a1))
	 wha2 <- which(bond.dat.single[,2] == as.character(a2))

	 wha1a2 <- intersect(wha1, wha2)

	 if(length(wha1a2) == 0)
		{wha1 <- which(bond.dat.single[,1] == as.character(a2))
	 	 wha2 <- which(bond.dat.single[,2] == as.character(a1))

	 	 wha1a2 <- intersect(wha1, wha2)
		}

	 bond.dat.single[wha1a2[1],4]
	}

detect_case <- function(h) # Determine the 'case' of the dihedral angle (see page 8899 of Mayo et al); 
	{a1 <- dihedral_array[h,2]; a2 <- dihedral_array[h,3] # Atoms in the bond
	 b.order <- get.dihedral.order(h)

	 t1 <- AtTypes[AtIndex[a1]]; t2 <- AtTypes[AtIndex[a2]] # Atom hybridizations

	 # Default is case A
	 case <- "A";

	 # Detect case B
	 if(grepl("_2", t1) & grepl("_3", t2) & b.order == "single" | grepl("_2", t2) & grepl("_3", t1) & b.order == "single")
		case <- "B"

	 if(grepl("_R", t1) & grepl("_3", t2) & b.order == "single" |	 grepl("_R", t2) & grepl("_3", t1) & b.order == "single")
		case <- "B"

	 # Detect case C
	 if(grepl("_2", t1) & grepl("_2", t2) & b.order == "double")
		case <- "C"

	 # Detect case D
	 if(grepl("_R", t1) & grepl("_R", t2)  & b.order == "resonance")
		case <- "D"

	 # Detect case E
	 
	 if(grepl("_2", t1) & grepl("_2", t2) & b.order == "single")
		case <- "E"

	 if(grepl("_R", t1) & grepl("_2", t2)  & b.order == "single" | grepl("_R", t2) & grepl("_2", t1)  & b.order == "single")
		case <- "E"
		
	 # Detect case F	
		
	 # Detect case H
	 
	 # Detect case I

	 if(t1 == "O_3" & t2 == "C_R" | t2 == "O_3" & t1 == "C_R")
		case <- "I"

	 if(t1 == "O_3" & t2 == "C_2" | t2 == "O_3" & t1 == "C_2")
		case <- "I"
		
	 if(t1 == "S_3" & t2 == "C_R" | t2 == "S_3" & t1 == "C_R")
		case <- "I"

	 if(t1 == "S_3" & t2 == "C_2" | t2 == "S_3" & t1 == "C_2")
		case <- "I"		
		
	 # Detect case J
	 
	 aI <- dihedral_array[h,1]; aK <- dihedral_array[h,4] # Atoms in the bond
	 
	 tI <- AtTypes[AtIndex[a1]]; tK <- AtTypes[AtIndex[a2]] # Atom hybridizations 
	 
	 if(grepl("_R", t1) & grepl("_3", t2) & !grepl("_2", tI) & !grepl("_R", tI) & b.order == "single" |	  
	    grepl("_R", t2) & grepl("_3", t1) & !grepl("_2", tK) & !grepl("_R", tK) & b.order == "single")
		case <- "J"	 
	 
	 if(grepl("_2", t1) & grepl("_3", t2) & !grepl("_2", tI) & !grepl("_R", tI) & b.order == "single" |	  
	    grepl("_2", t2) & grepl("_3", t1) & !grepl("_2", tK) & !grepl("_R", tK) & b.order == "single")
		case <- "J"	 	 
	 
	 # Output
	 
	 case
	}

if(length(wh.dihedral) > 0)
{
#dihedral_list <- sapply(wh.dihedral, dindex.get, USE.NAMES = FALSE)
	# Not giving the correct formatting!?

dihedral_list <- lapply(wh.dihedral, dindex.get)

dihedral_array <- dihedral_list[[1]] # This is the array of dihedral angles

if(length(dihedral_list) > 1)
{
for(j in 2:length(dihedral_list))
	dihedral_array <- rbind(dihedral_array, dihedral_list[[j]])
}

dihedral.dat <- array(0,c(nrow(dihedral_array), 4))
		# Column 1 = type, Column 2 = V_jk, Column 3 = n, Column = phi_jk (see page 8899 of Mayo et al).

for(k in 1:nrow(dihedral_array)) # Need to 'renormalize' to account for number possibilities; see pg. 8889 of Mayo et al.
	{# Determine number of dihedral bonds possible for central bond

	 a1 <- dihedral_array[k,2]
	 a2 <- dihedral_array[k,3]

	 wh1 <- which(dihedral_array[,2] == a1)
	 wh2 <- which(dihedral_array[,3] == a2) # Shouldn't need to account for reverse cases?

	 #wh3 <- which(dihedral_array[,2] == a2)
	 #wh4 <- which(dihedral_array[,3] == a1) # Shouldn't need to account for reverse cases?

	 n_renorm <- length(intersect(wh1,wh2))

	 casek <- detect_case(k)

	 if(casek == "A")
		{V_jk = 2.0/n_renorm; n_jk = 3; phi_jk = 180;
		}

	 if(casek == "B")
		{V_jk = 1.0/n_renorm; n_jk = 6; phi_jk = 0;
		}

	 if(casek == "C")
		{V_jk = 45/n_renorm; n_jk = 2; phi_jk = 180;
		}

	 if(casek == "D")
		{V_jk = 25/n_renorm; n_jk = 2; phi_jk = 180;
		}

	 if(casek == "E")
		{V_jk = 5/n_renorm; n_jk = 2; phi_jk = 180;
		}
		
	 if(casek == "F")
		{V_jk = 10/n_renorm; n_jk = 2; phi_jk = 180;
		}		
		
	 if(casek == "H")
		{V_jk = 2.0/n_renorm; n_jk = 2; phi_jk = 90;
		}

	 if(casek == "I")
		{V_jk = 2/n_renorm; n_jk = 2; phi_jk = 180;
		}
		
	 if(casek == "J")
		{V_jk = 2.0/n_renorm; n_jk = 3; phi_jk = 180;
		}

	 dihedral.dat[k,] <- c(casek, V_jk/2, n_jk, phi_jk)
	}

# Re-determinte the dihedral types, having accounted for renormalization

Dtypes_0 <- unique(dihedral.dat[,1])
dihedral.dat.1 <- dihedral.dat # This will be modified to account for the new types

for(k in 1:length(Dtypes_0))
	{whk <- which(dihedral.dat.1[,1] == Dtypes_0[k])
	 dihed.k <- dihedral.dat[whk,]
	 
	 if(length(whk) == 1)
		dihed.k <- t(dihed.k)
	 
	 Ktypes <- unique(dihed.k[,2])

	 # Make list of new types
	 NewTypes <- c()

	 for(h in 1:length(Ktypes))
		{NewTypes <- c(NewTypes, paste(Dtypes_0[k], h, sep = ""))

		 whh <- which(dihed.k[,2] == Ktypes[h])
		 dihed.k[whh,1] <- NewTypes[h]
		}

	 dihedral.dat.1[whk,] <- dihed.k
	}

dihedral.dat <- dihedral.dat.1 # For atom types

# Include additional molecules

if(nMol > 1)
	{dihedral.dat.base <- dihedral.dat
	 dihedral_array_base <- dihedral_array

	 for(j in 2:nMol)
		{dihedral_array_add <- dihedral_array_base + (j-1)*nAt
	 	 dihedral_array <- rbind(dihedral_array, dihedral_array_add)
		 dihedral.dat <- rbind(dihedral.dat, dihedral.dat.base)
		}
	}
}	

if(length(wh.dihedral) == 0)
	dihedral.dat <- array("",c(0,4))

## Determine impropers (for this molecule, only C_2 and C_R need to be kept planar)
## Is this correct? Check

whImp <- c(which(AtTypes == "C_2"), which(AtTypes == "C_R"), which(AtTypes == "N_2"))
AtImp <- AtIndex[whImp] # Atoms which carry an improper 

#AtImp <- c()

#for(j in 1:nMol)
#	{addj <- AtIndex[whImp] + (j-1) * nAt	
#	 AtImp <- c(AtImp, addj)
#	}

inv.dat <- array(0, c(length(AtImp), 5)) # First column is atom identity, second - fourth columns are identitues of three connected atoms, final column is force constant KL

whkeep <- c()

for(k in 1:length(AtImp))
	{wha1 <- which(bond.mat[,1] == AtImp[k])
	 neighbors <- bond.mat[wha1,2]

	 wha2 <- which(bond.mat[,2] == AtImp[k])
	 neighbors <- c(neighbors, bond.mat[wha2,1])

	 neighbors <- unique(neighbors)
	 
	 # Don't need - impropers should always involve 4 atoms!
	 #nRmn <- 3 - length(neighbors)
	 #neighbors <- c(neighbors, rep(0, nRmn))
	 
	 if(length(neighbors) == 3)
		{inv.dat[k,] <- c(AtImp[k], neighbors, 40)
		 whkeep <- c(whkeep, k)
		}
	}
	
inv.dat <- inv.dat[whkeep,]	

if(length(whkeep) == 1)
	inv.dat <- t(inv.dat)

if(nMol > 1)
	{inv.dat.base <- inv.dat

	 for(j in 2:nMol)
		{inv.dat.add <- inv.dat.base + (j-1)*nAt
		 inv.dat <- rbind(inv.dat, inv.dat.add)
		}
	}
	
############################################################################################
############### PART 2 - Create the LAMMPS .dat file #######################################
############################################################################################

## Set the bond and atom type data

AtomT <- AtTypesU # Atom types
BondT <- unique(bond.dat[,3]) # Bond types
AngleT <- unique(ang.dat[,1]) # Angle types
DihedralT <- unique(dihedral.dat[,1]) # Dihedral types
ImproperT <- 1 # Improper types

if(nWat > 0)
	{BondT <- unique(c(bond.dat[,3],bond.dat.h2o[,3]))
	 AngleT <- unique(c(ang.dat[,1],ang.dat.h2o[,1]))
	}
	
## Create the header of the .dat file

dat.header <- "LAMMPS data file"

dat.header <- c(dat.header, paste(" ", nMol*length(AtIndex) + 3*nWat, " atoms", sep = ""))
dat.header <- c(dat.header, paste(" ", nMol*nrow(bond.mat) + 2*nWat, " bonds", sep = ""))
dat.header <- c(dat.header, paste(" ", nrow(ang.dat) + nWat, " angles", sep = ""))
dat.header <- c(dat.header, paste(" ", nrow(dihedral.dat), " dihedrals", sep = ""))
dat.header <- c(dat.header, paste(" ", nrow(inv.dat), " impropers", sep = ""))

dat.header <- c(dat.header, paste(" ", length(AtomT), " atom types", sep = ""))
dat.header <- c(dat.header, paste(" ", length(BondT), " bond types", sep = ""))
dat.header <- c(dat.header, paste(" ", length(AngleT), " angle types", sep = ""))
dat.header <- c(dat.header, paste(" ", length(DihedralT), " dihedral types", sep = ""))
dat.header <- c(dat.header, paste(" ", length(ImproperT), " improper types", sep = ""))

dat.header <- c(dat.header, paste(c("", xlo, xhi, "xlo", "xhi"), collapse = " "))
dat.header <- c(dat.header, paste(c("", ylo, yhi, "ylo", "yhi"), collapse = " "))
dat.header <- c(dat.header, paste(c("", zlo, zhi, "zlo", "zhi"), collapse = " "))
dat.header <- c(dat.header, paste(c("", xy, xz, yz, "xy xz yz"), collapse = " "))

# Comments for atom types

dat.header <- c(dat.header, "")
dat.header <- c(dat.header, "# Atom types")
dat.header <- c(dat.header, "#")

count <- 1;
for(k in 1:length(AtomT))
	{addk <- paste("# ", count, " ", AtomT[k], sep = "")
	 dat.header <- c(dat.header, addk)

	 count <- count + 1
	}

# Comments for Bond Coeffs

dat.header <- c(dat.header, "")
dat.header <- c(dat.header, "# Bond Coeffs")
dat.header <- c(dat.header, "#")

count <- 1;
for(k in 1:length(BondT))
	{addk <- paste("# ", count, " ", BondT[k], sep = "")
	 dat.header <- c(dat.header, addk)

	 count <- count + 1
	}

# Comments for Angle Coeffs

dat.header <- c(dat.header, "")
dat.header <- c(dat.header, "# Angle Coeffs")
dat.header <- c(dat.header, "#")

count <- 1;
for(k in 1:length(AngleT))
	{addk <- paste("# ", count, " ", AngleT[k], sep = "")
	 dat.header <- c(dat.header, addk)

	 count <- count + 1
	}

# Comments for Dihedral Coeffs

dat.header <- c(dat.header, "")
dat.header <- c(dat.header, "# Dihedral Coeffs")
dat.header <- c(dat.header, "#")

if(length(DihedralT) > 0)
{
count <- 1;
for(k in 1:length(DihedralT))
	{addk <- paste("# ", count, " ", DihedralT[k], sep = "")
	 dat.header <- c(dat.header, addk)

	 count <- count + 1
	}
}

# Comments for Inversion Coeffs

dat.header <- c(dat.header, "")
dat.header <- c(dat.header, "# Inversion Coeffs")
dat.header <- c(dat.header, "#")

count <- 1;
for(k in 1:length(ImproperT))
	{addk <- paste("# ", count, " ", ImproperT[k], sep = "")
	 dat.header <- c(dat.header, addk)

	 count <- count + 1
	}

# Mass information

dat.header <- c(dat.header, "")
dat.header <- c(dat.header, " Masses")
dat.header <- c(dat.header, "")

for(k in 1:length(massT))
	{dat.header <- c(dat.header, paste(" ", k, " ", massT[k], sep = ""))
	}

# Atom information

dat.header <- c(dat.header, "")
dat.header <- c(dat.header, " Atoms # full")
dat.header <- c(dat.header, "")

for(j in 1:nMol) # For molecules
	{
	 for(k in 1:length(AtIndex))
		{tk <- AtTypes[k]

	 	 wht <- which(AtomT == tk)
	 	 whk <- which(Zcoord[,1] == k)

	 	 addk <- paste(c(k + (j-1)*nrow(Zcoord), j, wht, 0, ZcoordL[[j]][k,]), collapse = " ") # Zero is the charge; this may need to change ## CAREFuL! May need ZcoordL[[j]][whk,], if atom order not preserved
	 
	 	 dat.header <- c(dat.header, addk)
		}
	}

if(nWat > 0)
{
count <- nMol*nAt + 1 # For water

wt1 <- which(AtomT == "H__HB")
wt2 <- which(AtomT == "O_3")
wt3 <- wt1

for(j in 1:nWat)
	{addk1 <- paste(c(count, j + nMol, wt1, 0, WcoordL[[j]][1,]), collapse = " ")
	 addk2 <- paste(c(count + 1, j + nMol, wt2, 0, WcoordL[[j]][2,]), collapse = " ")
	 addk3 <- paste(c(count + 2, j + nMol, wt3, 0, WcoordL[[j]][3,]), collapse = " ")

	 dat.header <- c(dat.header, addk1, addk2, addk3)
	 
	 count <- count + 3
	}
}

# Bond information

dat.header <- c(dat.header, "")
dat.header <- c(dat.header, " Bonds")
dat.header <- c(dat.header, "")

if(nrow(bond.dat) > 0)
{
for(k in 1:nrow(bond.dat))
	{tk <- bond.dat[k,3]

	 wht <- which(BondT == tk)

	 addk <- paste(c(k, wht, as.numeric(bond.dat[k,1:2])), collapse = " ")

	 dat.header <- c(dat.header, addk)
	}
}

if(nWat > 0)
{
count <- nrow(bond.dat) + 1
addk2 <- which(BondT == "O_3*H__HB")

for(k in 1:nrow(bond.dat.h2o))
	{addk1 <- count;
	 addk3 <- bond.dat.h2o[k,1]
	 addk4 <- bond.dat.h2o[k,2]

	 addk <- paste(c(addk1, addk2, addk3, addk4), collapse = " ")

	 dat.header <- c(dat.header, addk)

	 count <- count + 1
	}
}

# Angle information

dat.header <- c(dat.header, "")
dat.header <- c(dat.header, " Angles")
dat.header <- c(dat.header, "")

if(nrow(ang.dat) > 0)
{
for(k in 1:nrow(ang.dat))
	{tk <- ang.dat[k,1]

	 wht <- which(AngleT == tk)

	 addk <- paste(c(k, wht, as.numeric(ang.dat[k,4:6])), collapse = " ")

	 dat.header <- c(dat.header, addk)
	}
}
	
if(nWat > 0)
{
count <- nrow(ang.dat) + 1 # For water molecules
addk2 <- which(AngleT == "O_3")

for(k in 1:nrow(ang.dat.h2o))
	{addk1 <- count;
	 addk3 <- ang.dat.h2o[k,4] 
	 addk4 <- ang.dat.h2o[k,5] 
	 addk5 <- ang.dat.h2o[k,6]

	 addk <- paste(c(addk1, addk2, addk3, addk4, addk5), collapse = " ")

	 dat.header <- c(dat.header, addk)

	 count <- count + 1
	}
}	

# Dihedral information

if(nrow(dihedral.dat) > 0)
{
dat.header <- c(dat.header, "")
dat.header <- c(dat.header, " Dihedrals")
dat.header <- c(dat.header, "")

for(k in 1:nrow(dihedral.dat))
	{tk <- dihedral.dat[k,1]

	 wht <- which(DihedralT == tk)

	 addk <- paste(c(k, wht, dihedral_array[k,]), collapse = " ")

	 dat.header <- c(dat.header, addk)
	}
}

# Improper information

if(nrow(inv.dat) > 0)
{
dat.header <- c(dat.header, "")
dat.header <- c(dat.header, " Impropers")
dat.header <- c(dat.header, "")

for(k in 1:nrow(inv.dat))
	{tk <- 1

	 wht <- which(ImproperT == tk)

	 addk <- paste(c(k, wht, inv.dat[k,1:4]), collapse = " ")

	 dat.header <- c(dat.header, addk)
	}
}

## Write the data file

if(autowrite.dat)
	writeLines(dat.header, dat.output)

############################################################################################
############### PART 3 - Create the LAMMPS .in file ########################################
############################################################################################

in.header <- "# Test LAMMPS input file for single DeFDe"
in.header <- c(in.header, "")
in.header <- c(in.header, "# Initialize simulation")

## Basic simulation settings

in.header <- c(in.header, "")

in.header <- c(in.header, "clear")
in.header <- c(in.header, "")

in.header <- c(in.header, "units real")
in.header <- c(in.header, "atom_style full")
in.header <- c(in.header, "boundary p p p")
in.header <- c(in.header, "dielectric 1") # Just copying the example
in.header <- c(in.header, "special_bonds dreiding") 
in.header <- c(in.header, "")

in.header <- c(in.header, "bond_style harmonic")
in.header <- c(in.header, "angle_style harmonic")
in.header <- c(in.header, "dihedral_style harmonic")
in.header <- c(in.header, "improper_style umbrella")

# pair_style for hydrogen bonding
in.header <- c(in.header, "pair_style hybrid/overlay hbond/dreiding/lj 4 9.0 11.0 90 lj/cut 8.50000") # lj/cut should be OK for when charges are absent. 9 and 11 refer to the inner and outer cut-offs for the hydrogen bonding. See https://lammps.sandia.gov/doc/pair_hbond_dreiding.html

# pair_style for non-hydrogen bonding
#in.header <- c(in.header, "pair_style lj/cut 8.50000")

## Read data file

in.header <- c(in.header, "")
in.header <- c(in.header, "# Read data file")
in.header <- c(in.header, paste("read_data ", dat.output, sep = ""))

## DREIDING potential information

in.header <- c(in.header, "")
in.header <- c(in.header, "# Dreiding potential information")

# Bond information

in.header <- c(in.header, "")

for(k in 1:length(BondT))
	{tyk <- BondT[k]

	 bond.dat.k <- which(bond.dat[,3] == tyk)[1]
	 param.k <- as.numeric(bond.dat[bond.dat.k,5:6])
	
	 if(is.na(param.k[1])) # Water bond
		param.k <- as.numeric(bond.dat.h2o[1,5:6])
	
	 addk <- paste(c("bond_coeff", k, 0.5*param.k[2], param.k[1]), collapse = " ")

	 in.header <- c(in.header, addk)	 
	}

# Angle information

in.header <- c(in.header, "")

for(k in 1:length(AngleT))
	{tyk <- AngleT[k]

	 ang.dat.k <- which(ang.dat[,1] == tyk)[1]
	 param.k <- as.numeric(ang.dat[ang.dat.k,2:3])
	 
	 if(is.na(param.k[1]))
		param.k <- as.numeric(ang.dat.h2o[1,2:3])

	 addk <- paste(c("angle_coeff", k, 0.5*param.k[1], param.k[2]), collapse = " ")

	 in.header <- c(in.header, addk)
	}

# Dihedral information

in.header <- c(in.header, "")

if(length(DihedralT) > 0)
{
for(k in 1:length(DihedralT))
	{tyk <- DihedralT[k]

	 dih.dat.k <- which(dihedral.dat[,1] == tyk)[1]
	 param.k <- c(as.numeric(dihedral.dat[dih.dat.k,2]), as.numeric(dihedral.dat[dih.dat.k,3])) # Is the *0.5 necessary? No.
	 dadd <- -1

	 #if(strsplit(tyk, "")[[1]][1] == "B")
	 #	dadd <- -1 # Case H will require special treatment again.

	 addk <- paste(c("dihedral_coeff", k, param.k[1], dadd, param.k[2]), collapse = " ")

	 in.header <- c(in.header, addk)
	}
}

# Inproper information (only one type of improper is assumed at the moment)

in.header <- c(in.header, "")

addk <- paste(c("improper_coeff", 1, 40, 0), collapse = " ")

in.header <- c(in.header, addk)

# Pairwise interaction data

get.vdw.param <- function(i,j) # Get the vdW parameter for atom interaction i and j
	{Di <- as.numeric(vdW.dat[i,3])
	 Dj <- as.numeric(vdW.dat[j,3])

	 ri <- as.numeric(vdW.dat[i,2])
	 rj <- as.numeric(vdW.dat[j,2])

	 # Combine parameters as geometric means
	 Dij <- sqrt(Di*Dj)
	 rij <- sqrt(ri*rj)

	 c(Dij, rij)
	}

in.header <- c(in.header, "")

for(i in 1:length(AtomT))
	for(j in i:length(AtomT))
		{tyi <- AtomT[i]
		 tyj <- AtomT[j]

		 whi <- which(vdW.dat[,1] == tyi)
		 whj <- which(vdW.dat[,1] == tyj)

		 DijRij <- get.vdw.param(whi,whj)

		 Dij <- DijRij[1]
		 Rij <- DijRij[2] # These look rather long... Check definition

		 addk <- paste(c("pair_coeff", i, j, "lj/cut", Dij, Rij), collapse = " ") # When hydrogen bonding is present, lj/cut should be written explicit
		 
		 #addk <- paste(c("pair_coeff", i, j, Dij, Rij), collapse = " ") # When hydrogen bonding is not present

		 in.header <- c(in.header, addk)
		}

# Add the H-bond line (later)

whN_2 <- which(grepl("N_2", AtomT)) # Obtain the sp2 nitrogen atoms
whO_3 <- which(grepl("O_3", AtomT)) # Obtain the sp3 oxygen atoms
whO_2 <- which(grepl("O_2", AtomT)) # Obtain the sp2 oxygen atoms
whF <- which(grepl("F_", AtomT)) # Obtain the fluorine atoms
whHB <- which(AtomT == "H__HB") # Hydrogen bonding atom

# Only the O_3 atom on water is the donor 

whDonor <- sort(unique(c(whN_2, whO_3)))
whAcceptor <- sort(unique(c(whN_2, whO_2, whO_3, whF)))

# Sort the acceptors-donors, so that small appears before big

whE <- expand.grid(whAcceptor, whDonor)

da_list <- as.matrix(whE);
rownames(da_list) <- NULL
colnames(da_list) <- NULL

for(k in 1:nrow(da_list))
	{wh1 <- da_list[k,1]; wh2 <- da_list[k,2]
	 
	 don = "j"
	 if(wh1 > wh2)
	 	{a <- wh1
		 wh1 <- wh2
		 wh2 <- a
		 don = "i"
		}

	 addk1 <- paste("pair_coeff ", wh1, " ", wh2, " hbond/dreiding/lj ", whHB, " ", don, " 9.50 2.75 4 9.0 11.0 90", sep = "")
	 
	 in.header <- c(in.header, addk1)
	}

# Calculation parameters

# Calculation execution

in.header <- c(in.header,
"",
"neighbor 2 bin",
"neigh_modify every 2 one 10000",
"",
"# Gemoetric relaxation",
"fix 1 all box/relax iso 1.0 vmax 0.001",
"minimize 1.0e-4 1.0e-6 10000 10000",
"",
"# MD simulation",
"velocity all create 300.0 12345 rot yes dist gaussian",
"timestep 1.0",
"fix 2 all nvt temp 300.0 300.0 1",
"",
"# Compute interaction energy",
paste("group WAT molecule <>", length(ZcoordL) + 1, length(ZcoordL) + length(WcoordL), sep = " "),
"group MOL molecule 1",
"",
"compute w WAT group/group MOL pair yes",
"fix 3 all ave/time 10 5 100 c_w file outa",
"",
"# Outputs",
"thermo 100 # Output thermodynamic data every nth step",
"#dump mf1 all xyz 10000 mol_*.xyz ",
"",
"# Run",
"run 1000000",
"",
"# Output the system volume",
"variable v equal vol",
"print 'The system volume is $v'",
"",
"# Outputs",
"variable A equal xhi-xlo",
"variable B equal yhi-ylo",
"variable C equal zhi-zlo",
"variable D equal xy",
"variable E equal xz",
"variable F equal yz",
"print 'Lattice vectors'",
"print '$A 0 0'",
"print '$D $B 0'",
"print '$E $F $C'"
)

## Write the data file

if(autowrite.in)
	writeLines(in.header, in.output)







