#!/usr/bin/env perl

use strict;
use warnings;
use lib '/data/store/spb41/apps/perl/lib64/perl5';
use constant PI=>3.1415926;
use POSIX qw(ceil floor);
use File::Spec;
use Getopt::Long;
# We want the following to be specified on the command line, with 
# a default.
# my $Directory="/local1/spb41/test/";
my $Directory="";
my $GlassFile="";
my $GlassPart='';
#"/home/spb41/data2/small-sim/";
#Simulation parameters
my $seed='';
#my $seed=128510;
my $boxsize='';
my $npart='';
my $kspace='';
# 1 is Eisenstein & Hu
# 2 is some tabulated power spectrum written by Volker
# 3 is a CAMB transfer function, with different transfer functions for different species, but P(k) = A k^n_s
# 4 is a spline, CAMB transfer function, different transfers for different species
# 5 is a spline, CAMB transfer function, *total* transfer for all species.
my $Spectrum=4;
#Cosmological parameters
my $Omega_M='';
my $Omega_L='';
my $Omega_B='';
my $M_Nu='';
my $Redshift='';
my $NS='';
my $A_prim='';
my $hub='';
my $Sigma_8=0.8;
my $Prefix='';
#my $Sigma_8=2.4876783E-02;#Should be sigma 8 given by CAMB. Isn't actually used. 

#For printing help
my $help='';
my @options=@ARGV;
GetOptions ('seed:i'=>\$seed, 'box:f'=>\$boxsize, 'npart:i'=>\$npart,'om:f'=>\$Omega_M, 'ol:f'=>\$Omega_L,
        'ob:f'=>\$Omega_B,'redshift:f'=>\$Redshift, 'output=s'=>\$Directory, 
        'glass=s'=>\$GlassFile, 'nglass:i'=>\$GlassPart,'mnu:f'=>\$M_Nu, 'help'=>\$help,
         'hub:f'=>\$hub, "prefix:s"=>\$Prefix, 'kspace'=>\$kspace, 'ns:f'=>\$NS,
         'as:f'=>\$A_prim
                ) or die "Failed Options";
if($help){
        print "Usage: 
               --seed 250       => Seed value for the noise.
               --box  60        => Boxsize in Mpc/h
               --npart 400      => Number of particles.
               --om   f         => Omega_matter
               --ol   f         => Omega_Lambda
               --ob   f         => Omega_Baryons
               --mnu  f         => M_nu (single flavour)
               --prefix str     => Prefix for IC file
               --redshift n     => Starting redshift
               --output  path   => Directory to output files to
               --glass   path   => Glass file to use
               --nglass  100    => Number of particles in glass file
               --help           => This message";
               die;
       }
       #Parameters are WMAP 7-year.
if(!$Omega_M){$Omega_M=0.222;}
if(!$Omega_B){$Omega_B=0.0449;}
if(!$Omega_L){$Omega_L = 1-$Omega_M-$Omega_B;}
if(!$hub){$hub=0.70;}
if(!$Redshift){$Redshift=99;}
if(!$M_Nu){$M_Nu = 0;}
if(!$kspace){$kspace = 0;}
if(!$NS){$NS=1.0;}
if(!$A_prim){$A_prim=2.43e-9;}
my $Omega_Nu = 3*$M_Nu / 93.14/$hub/$hub;
$Omega_M -= $Omega_Nu;
#Keep the seed the same as in the best-fit case.
if(!$seed){$seed=181170;}
if(!$boxsize){$boxsize=60;}
if(!$GlassPart){$GlassPart=100;}
if(!$npart){$npart=400;}
if(!$Prefix){$Prefix="ics";}
#Add the prefix onto the end of the output directory. 
$Directory= "$Directory";

#Paths: Try to find them with globbing.
my $NGenIC="$ENV{HOME}/N-GenIC/N-GenIC";
unless (-e $NGenIC){
        my @lst = glob("$ENV{HOME}/*/N-GenIC/N-GenIC");
        $NGenIC = $lst[0];
}
my $CAMB="$ENV{HOME}/cosmomc/camb/camb";
unless (-e $CAMB){
        my @lst = glob("$ENV{HOME}/*/cosmomc/camb/camb");
        if( $#lst > -1){
            $CAMB = $lst[0];
        }
        else{
            @lst = glob("$ENV{HOME}/*/*/cosmomc/camb/camb");
            ( $#lst > -1) or die "Could not find CAMB!\n";
            $CAMB = $lst[0];
        }
}
my $class="$ENV{HOME}/codes/class/class";
unless (-e $class){
        my @lst = glob("$ENV{HOME}/*/class/class");
        $class = $lst[0];
}
my $genpk= "$ENV{HOME}/C-post-process/gen-pk";
unless (-e $genpk){
        my @lst = glob("$ENV{HOME}/*/C-post-process/gen-pk");
        $genpk = $lst[0];
}
#Files
#If we have more than 4GB of particle data, 
#the size of the POS area will overflow.
#Also we don't want N-GenICs to take too long.
my $NumFiles=ceil($npart**3*3/2**29);
if($Omega_B > 0 || $M_Nu > 0 ){
        $NumFiles=ceil($npart**3*3/2**28);
}
# We don't ever want to use an odd number of files or processors.
if($NumFiles > 1){
        $NumFiles = 2*ceil($NumFiles/2);
}
if($NumFiles < 4){
        $NumFiles =4;
}
my $GlassTileFac=$npart/$GlassPart;
my $ICFile=$Prefix."_$npart-$boxsize-z$Redshift-nu$M_Nu.dat";
my $NGenParams="_spb41-gen.param";
my $NGenDefParams=$NGenIC;
$NGenDefParams =~ s!/[\w-]*$!/ics_nu_default.param!;
my $CAMBParams="_spb41-camb-params.ini";
my $CAMBDefParams = $CAMB;
$CAMBDefParams = "/home/spb/codes/numatpowscripts/paramfiles/default-params.ini"; #=~ s!/[\w-]*$!/paramfiles/default-params.ini!;
my $classprecision =$class;
$classprecision =~ s/class$/pk_ref.pre/;
my $classdefparams = $class;
$classdefparams =~ s/class$/explanatory.ini/;
my $Pkestimate="pk-init-$npart-$boxsize-z$Redshift-$seed";
#my $TransferFile=$Prefix."_transfer_$Redshift.dat";
my $TransferFile=$Prefix."_z1_tk.dat";
my $PYPlotScript="_plot-init.py";

# $CAMB="true";
# $NGenIC="true";
# $genpk= "true";
#Output the parameter file for CAMB.

#Output the initial conditions file for N-GenICs.
$CAMBParams=$Directory."/".$CAMBParams;
#Generate later redshifts
my @Redshifts = ($Redshift);
my @rrr = (49,24,9,4,3,2,1,0);
foreach(@rrr)
{
        if($_ < $Redshift){
                push @Redshifts,$_;
        }
}
# paramfile newparams output_root omega_nu omega_b omega_cdm hubble redshift
#gen_camb_file($CAMBDefParams,$CAMBParams,$Directory."/".$Prefix,$Omega_Nu,$Omega_B, $Omega_M, $hub, $NS, $A_prim, $Redshift);
gen_class_file($classdefparams,$CAMBParams,$Directory."/".$Prefix,$M_Nu,$Omega_B, $Omega_M, $hub, $NS, $A_prim, @Redshifts);
print "Running CLASS...\n";
#print_run("$CAMB $CAMBParams");
print_run("$class $CAMBParams $classprecision");
print "Done Running CLASS, running N-GenICs.\n";
$NGenParams=$Directory."/".$NGenParams;
# Output the initial conditions file for N-GenICs.
# paramfile newparams directory transferfile icfile glassfile glasspart npart box omega_nu omega_b omega_m hubble redshift kspace
gen_ngen_file($NGenDefParams , $NGenParams,$Directory, $TransferFile,$ICFile,$GlassFile, $GlassPart,$npart, $boxsize,$Omega_Nu, $Omega_B, $Omega_M, $hub, $Redshift, $kspace, $NumFiles);
print_run("$NGenIC $NGenParams");

#If we are using Kspace neutrinos, we also want to generate CAMB tables for all redshifts along the way.
#Not for the new style kspace neutrinos
if( 0  and $kspace){
        my @Redshifts;
        for(my $i = log(1+$Redshift); $i > 0; $i-=1./50){
                push @Redshifts, exp($i)-1;
        }
        push @Redshifts,0;
        print "Generating CAMB_TABLES\n";
        `mkdir -p $Directory/CAMB_TABLES/`;
        my $IntParams = $Directory."/CAMB_TABLES/_int-camb-params.ini";
        gen_camb_file($CAMBDefParams,$IntParams,$Directory."/CAMB_TABLES/tab",$Omega_Nu,$Omega_B, $Omega_M, $hub, $NS, $A_prim, @Redshifts);
        print_run("$CAMB $IntParams");
}

print "Done generating initial conditions, now finding P(k).\n";

print_run("$genpk -i $Directory/$ICFile -o $Directory");
#Output file describing run.
my $txtfile =$Directory."/".$Prefix.".txt";
open(my $outhandle, ">", $txtfile) or 
      die "Can't open file $txtfile for writing!";
print $outhandle "
perl $0 @options
";
close($outhandle);

# pyscript datafile dir O_M box hub redshift
gen_plot_script($PYPlotScript, $ICFile, $Directory,$Prefix,$Omega_M, $Omega_Nu, $boxsize, $hub, $Redshift,$Omega_B, $kspace);
#Execute the script
print `python $PYPlotScript`;

# paramfile newparams output_root omega_nu omega_b omega_cdm hubble redshift
sub gen_camb_file{
        my $paramfile=shift;
        my $newparams=shift;
        my $root = shift;
        my $o_nu=shift; 
        my $o_b = shift;
        my $o_cdm=shift;
        if($o_b < 0.001){ 
                $o_b = 0.05;
                $o_cdm -=$o_b;
        }
        my $hub = 100*(shift);
        my $NS = shift;
        my $A_prim = shift;
        my @red = @_;
        my $mass_nu;
        my $nomass_nu;
        if($o_nu == 0){
                $mass_nu = 0;
                $nomass_nu = 3.04;
        }else{
                $mass_nu = 3.04;
                $nomass_nu = 0;
        }
        #Read in template parameter file
        open(my $INHAND, "<","$paramfile") or die "Could not open $paramfile for reading!";
        open(my $OUTHAND, ">","$newparams") or die "Could not open $newparams for writing!";
        while(<$INHAND>){
                s/^\s*output_root\s*=\s*[\w\/\.-]*/output_root= $root/i;
                #Set neutrino mass
                s/^\s*use_physical\s*=\s*[\w\/.-]*/use_physical = F/i;
                s/^\s*omega_cdm\s*=\s*[\w\/.-]*/omega_cdm=$o_cdm/i;
                s/^\s*omega_baryon\s*=\s*[\w\/.-]*/omega_baryon=$o_b/i;
                s/^\s*omega_lambda\s*=\s*[\w\/.-]*/omega_lambda=$Omega_L/i;
                s/^\s*omega_neutrino\s*=\s*[\w\/.-]*/omega_neutrino=$o_nu/i;
                s/^\s*hubble\s*=\s*[\w\/.-]*/hubble = $hub/i;
                #Set things we always need here
                s/^\s*massless_neutrinos\s*=\s*[\w\/.-]*/massless_neutrinos = $nomass_nu/i;
                s/^\s*massive_neutrinos\s*=\s*[\w\/.-]*/massive_neutrinos = $mass_nu/i;
                s/^\s*get_transfer\s*=\s*[\w\/.-]*/get_transfer = T/i;
                s/^\s*do_nonlinear\s*=\s*[\w\/.-]*/do_nonlinear = 0/i;
                #Set initial conditions; scalar_amp is to give sigma_8 = 0.878 at z=0 with nu=0.
                #Pivot irrelevant as n_s = 1
                s/^\s*initial_power_num\s*=\s*[\w\/.-]*/initial_power_num = 1/i;
                s/^\s*scalar_amp\(1\)\s*=\s*[\w\/.-]*/scalar_amp(1) = $A_prim/i;
                s/^\s*scalar_spectral_index\(1\)\s*=\s*[\w\/.-]*/scalar_spectral_index(1) = $NS/i;
                s/^\s*scalar_nrun\(1\)\s*=\s*[\w\/.-]*/scalar_nrun(1) = 0/i;
		#Set the pivot scale for scalar and tensor modes to be WMAP value
                s/^\s*pivot_scalar\s*=\s*[\w\/.-]*/pivot_scalar = 2e-3/i;
                s/^\s*pivot_tensor\s*=\s*[\w\/.-]*/pivot_tensor = 2e-3/i;
                #Set up output
                s/^\s*transfer_kmax\s*=\s*[\w\/.-]*/transfer_kmax = 200/i;
                s/^\s*transfer_k_per_logint \s*=\s*[\w\/.-]*/transfer_k_per_logint  = 30/i;
                my $nout=$#red+1;
                s/^\s*transfer_num_redshifts\s*=\s*[\w\/.-]*/transfer_num_redshifts = $nout/i;
                #It is EXTREMELY IMPORTANT that this is set to F, because the 
                #calculations done in N-GenICs and Gadget assume that a matterpow 
                #row has the same k as the equivalent row in the transfer function!
                s/^\s*transfer_interp_matterpower\s*=\s*[\w\/.-]*/transfer_interp_matterpower = F/i;
                #Output files set later.
                s/^\s*transfer_redshift\([0-9]\)\s*=\s*[\w\/.-]*//i;
                s/^\s*transfer_filename\([0-9]\)\s*=\s*[\w\/.-]*//i;
                s/^\s*transfer_matterpower\([0-9]\)\s*=\s*[\w\/.-]*//i;
                #Write to new file.
                print $OUTHAND $_;
        }
        #Set output files
        print $OUTHAND "\n#Transfer output files\n";
        for(my $i=0; $i<=$#red; $i++){
                print $OUTHAND "transfer_redshift(".($i+1).") = $red[$i]\n";
                print $OUTHAND "transfer_filename(".($i+1).") = transfer_$red[$i].dat\n";
                print $OUTHAND "transfer_matterpower(".($i+1).") = matterpow_$red[$i].dat\n";
        }
        close($INHAND);
        close($OUTHAND);
}

# paramfile newparams output_root omega_nu omega_b omega_cdm hubble redshift
sub gen_class_file{
        my $paramfile=shift;
        my $newparams=shift;
        my $root = shift;
        $root = $root."_";
        my $M_Nu=shift; 
        my $o_b = shift;
        my $o_cdm=shift;
        if($o_b < 0.001){ 
                $o_b = 0.05;
                $o_cdm -=$o_b;
        }
        my $hub = shift;
        my $NS = shift;
        my $A_prim = shift;
        my @red = @_;
        my $mass_nu;
        my $nomass_nu;
        if($M_Nu == 0){
                $mass_nu = 0;
                $nomass_nu = 3.04;
        }else{
                $mass_nu = 3;
                $nomass_nu = 0;
        }
        my $redstr=join(", ",@red);
        #Read in template parameter file
        open(my $INHAND, "<","$paramfile") or die "Could not open $paramfile for reading!";
        open(my $OUTHAND, ">","$newparams") or die "Could not open $newparams for writing!";
        while(<$INHAND>){
                s/^\s*root\s*=\s*[\w\/\.-]*/root = $root/i;
                #Set neutrino mass
                s/^\s*Omega_cdm\s*=\s*[\w\/.-]*/Omega_cdm=$o_cdm/i;
                s/^\s*Omega_b\s*=\s*[\w\/.-]*/Omega_b=$o_b/i;
                s/^\s*Omega_Lambda\s*=\s*[\w\/.-]*/Omega_Lambda=$Omega_L/i;
                s/^\s*m_ncdm\s*=\s*[\w\/., -]*/m_ncdm=$M_Nu,$M_Nu,$M_Nu/i;
                #Class does not support non-flat models
                s/^\s*Omega_k\s*=\s*[\w\/.,-]*/Omega_k=0/i;
                s/^\s*h\s*=\s*[\w\/.-]*/h = $hub/i;
                #Set things we always need here
                s/^\s*N_eff\s*=\s*[\w\/.-]*/N_eff = $nomass_nu/i;
                s/^\s*N_ncdm\s*=\s*[\w\/.-]*/N_ncdm = $mass_nu/i;
                s/^\s*output\s*=\s*[\w\/.,-]*/output = mPk,mTk/i;
                s/^\s*non linear\s*=\s*[\w\/.-]*/non linear = no/i;
                s/^\s*format\s*=\s*[\w\/.-]*/format = camb/i;
                s/^\s*headers\s*=\s*[\w\/.-]*/headers= no/i;
                s/^\s*YHe\s*=\s*[\w\/.-]*/YHe = 0.24/i;
                #Set initial conditions; scalar_amp is to give sigma_8 = 0.878 at z=0 with nu=0.
                s/^\sA_s\s*=\s*[\w\/.-]*/A_s = $A_prim/i;
                s/^\s*n_s\s*=\s*[\w\/.-]*/n_s = $NS/i;
                s/^\s*alpha_s\s*=\s*[\w\/.-]*/alpha_s = 0./i;
		#Set the pivot scale for scalar and tensor modes to be WMAP value
                s/^\s*k_pivot\s*=\s*[\w\/.-]*/k_pivot = 2e-3/i;
                #Set up output
                s!^\s*P_k_max_h/Mpc\s*=\s*[\w\/.-]*!P_k_max_h/Mpc = 100!i;
                s/^\s*transfer_k_per_logint \s*=\s*[\w\/.-]*/transfer_k_per_logint  = 30/i;
                #It is EXTREMELY IMPORTANT that this is set to F, because the 
                #calculations done in N-GenICs and Gadget assume that a matterpow 
                #row has the same k as the equivalent row in the transfer function!
                #s/^\s*transfer_interp_matterpower\s*=\s*[\w\/.-]*/transfer_interp_matterpower = F/i;
                #Output files set later.
                s/^\s*z_pk\s*=\s*[\w\/.-]*/z_pk = $redstr/i;
                #Write to new file.
                print $OUTHAND $_;
        }
        close($INHAND);
        close($OUTHAND);
}

# paramfile newparams directory transferfile icfile glassfile glasspart npart box omega_nu omega_b omega_m hubble redshift kspace
sub gen_ngen_file{
        my $NGenDefParams = shift;
        my $NGenParams=shift;
        my $Directory = shift;
        my $TransferFile=shift;
        my $MatterFile = $TransferFile;
        #$MatterFile=~ s/_transfer_/_matterpow_/;
        $MatterFile=~ s/tk/pk/;
        my $ICFile = shift;
        my $GlassFile= shift;
        my $GlassPart = shift;
        my $npart=shift;
        my $nmesh=3*$npart/2;
        my $GlassTileFac=$npart/$GlassPart;
        my $box=1000*shift; 
        my $O_nu = shift;
        my $Omega_B = shift;
        my $Omega_M=(shift) + $Omega_B +$O_nu;
        my $hub = shift;
        my $red = shift;
        my $kspace = shift;
        my $NumFiles=shift;
        my $nu = $O_nu*1.0 > 0 ? 1 : 0;
        #Subtract massless neutrinos
        if (!$nu){
                $Omega_M -=3.38663e-05;
        }
        my $nu_mass = 93.14*$hub*$hub*$O_nu;
        open(my $INHAND, "<","$NGenDefParams") or die "Could not open $NGenDefParams for reading!";
        open(my $OUTHAND, ">","$NGenParams") or die "Could not open $NGenParams for writing!";
        while(<$INHAND>){
                s/^\s*Nsample\s+[0-9]*/Nsample $npart/i;
                s/^\s*Nmesh\s+[0-9]*/Nmesh $nmesh/i;
                s!^\s*OutputDir\s+[/\w\.=-]*! OutputDir  $Directory!i;
                s!^\s*FileWithTransfer\s+[/\.\w=]*! FileWithTransfer  $Directory/$TransferFile!i;
                s!^\s*FileWithInputSpectrum\s+[/\.\w=]*! FileWithInputSpectrum  $Directory/$MatterFile!i;
                s!^\s*FileBase\s+[/\w]*! FileBase  $ICFile!i;
                s!^\s*GlassFile\s+[/\w\.=]*! GlassFile  $GlassFile!i;
                s!^\s*GlassTileFac\s+[/\w]*! GlassTileFac  $GlassTileFac!i;
                s!^\s*Box\s+[/\w]*! Box  $box!i;
                s!^\s*Seed\s+[/\w]*! Seed  $seed!i;
                s!^\s*Omega\s+[/\w\.]*! Omega  $Omega_M!i;
                s!^\s*OmegaDM_2ndSpecies\s+[/\w\.]*! OmegaDM_2ndSpecies  $O_nu!i;
                s!^\s*OmegaBaryon\s+[/\w\.]*! OmegaBaryon  $Omega_B!i;
                s!^\s*WhichSpectrum\s+[/\w]*! WhichSpectrum 2!i;
                s!^\s*OmegaLambda\s+[/\w\.]*! OmegaLambda  $Omega_L!i;
                s!^\s*NumFiles\s+[/\w\.]*! NumFiles  $NumFiles!i;
                s!^\s*Redshift\s+[/\w\.]*! Redshift  $red!i;
                s!^\s*HubbleParam\s+[/\w\.]*! HubbleParam  $hub!i;
                s!^\s*NU_On\s+[/\w]*! NU_On  $nu!i;
                #Do not renormalise the input spectrum
                s!^\s*ReNormalizeInputSpectrum\s+[/\w]*!ReNormalizeInputSpectrum 0!i;
                s!^\s*NU_KSPACE\s+[/\w]*! NU_KSPACE  $kspace!i;
                s!^\s*NU_Vtherm_On\s+[/\w]*! NU_Vtherm_On  1!i;
                s!^\s*NU_PartMass_in_ev\s+[/\w\.]*! NU_PartMass_in_ev  $nu_mass!i;
                print $OUTHAND $_;
        }
        close($INHAND);
        close($OUTHAND);
}

sub print_run{
        my $cmd=shift;
        $cmd=$cmd." |";
        my $PIPEHANDLE;
        open($PIPEHANDLE, $cmd);
        while(<$PIPEHANDLE>){
                print $_;
        }
        close($PIPEHANDLE);
        if($?){
           die "Problem running $cmd\n";
           }
}

# pyscript datafile dir O_M box hub redshift
sub gen_plot_script{
#Now we want to generate a script and plot all of these things. 
my $PyPlotScript= shift;
my $DataFile = shift;
my $Directory= shift;
my $Prefix = shift;
my $Omega_M = shift;
my $Omega_Nu = shift;
my $boxsize = shift;
my $hub = shift;
my $Redshift = shift;
my $Omega_B = shift;
my $kspace = shift;
$PYPlotScript =$Directory."/".$PYPlotScript;
open($outhandle, ">", $PYPlotScript) or 
      die "Can't open parameter file $PYPlotScript for writing!";
(my $la, my $la2,my $f) = File::Spec->splitpath($Directory."/".$DataFile);
print $outhandle 
'#!/usr/bin/python2

"""
Plot P(k)
"""

import numpy
import scipy.interpolate
import math
import matplotlib
matplotlib.use(\'agg\')
from matplotlib.pyplot import *


def plot_power(box,filename_dm,filename_b,camb_filename,redshift,hub,omegab,omegam):
        scale=2.0*math.pi/box
        camb_transfer=numpy.loadtxt(camb_filename)
        sigma=2.0
        #Adjust Fourier convention.
        k=camb_transfer[:,0]*hub
        #NOW THERE IS NO h in the T anywhere.
        Pk=camb_transfer[:,1]
        #^2*2*!PI^2*2.4e-9*k*hub^3
        ylabel("$P(k) /(Mpc^3 h^{-3})$")
        xlabel("$k /(h MPc^{-1})$")
        title("Power spectrum at z="+str(redshift))
        xlim(0.01, 100)
        loglog(k/hub, Pk, linestyle="--")
        pkc=numpy.loadtxt(filename_dm)
#        final=numpy.where(pkc[:,0] < pkc[-1,0]/1.0)
#       pkc=pkc[final]
        k_gadget=(pkc[1:,0])*scale
        # The factor of 2\pi^3 corrects for the units difference between gadget and camb
        Pk_gadget_dm=(2*math.pi)**3*pkc[1:,1]/scale**3
        if ( omegab > 0 ):
                pk_b=numpy.loadtxt(filename_b)
                #               pk_b=pk_b[final]
                Pk_gadget_b=(2*math.pi)**3*pk_b[1:,1]/scale**3
        else:
                Pk_gadget_b =0

        Pk_gadget=((omegam)*Pk_gadget_dm+omegab*Pk_gadget_b)/(omegam+omegab)
        samp_err=pkc[1:,2]
        sqrt_err=numpy.array(map(math.sqrt, samp_err))
        loglog(k_gadget,Pk_gadget, color="black", linewidth="1.5")
        loglog(k_gadget,Pk_gadget*(1+sigma*(2.0/sqrt_err+1.0/samp_err)),linestyle="-.",color="black")
        loglog(k_gadget,Pk_gadget*(1-sigma*(2.0/sqrt_err+1.0/samp_err)),linestyle="-.",color="black")
        xlim(0.01,k_gadget[-1]*1.1)
        return (k_gadget,Pk_gadget)

';
if(!$kspace and $Omega_Nu > 0){
print $outhandle
"
plot_power($boxsize,'$Directory/PK-DM-".$f."','$Directory/PK-nu-".$f."','$Directory/$Prefix"."_z1_pk.dat','$Redshift',$hub,$Omega_Nu,$Omega_M)
savefig('$Directory/$Pkestimate.png')
clf()
plot_power($boxsize,'$Directory/PK-DM-".$f."','$Directory/PK-nu-".$f."','$Directory/$Prefix"."_z1_pk.dat','$Redshift',$hub,$Omega_Nu,0)
savefig('$Directory/$Pkestimate-bar.png')
clf()";
}

if($Omega_B > 0){
print $outhandle
"
plot_power($boxsize,'$Directory/PK-DM-".$f."','$Directory/PK-by-".$f."','$Directory/$Prefix"."_z1_pk.dat','$Redshift',$hub,$Omega_B,$Omega_M)
savefig('$Directory/$Pkestimate.png')
clf()
plot_power($boxsize,'$Directory/PK-DM-".$f."','$Directory/PK-by-".$f."','$Directory/$Prefix"."_z1_pk.dat','$Redshift',$hub,$Omega_B,0)
savefig('$Directory/$Pkestimate-bar.png')
clf()";
}

print $outhandle
"
plot_power($boxsize,'$Directory/PK-DM-".$f."','$Directory/PK-DM-".$f."','$Directory/$Prefix"."_z1_pk.dat','$Redshift',$hub,0,$Omega_M)
savefig('$Directory/$Pkestimate-DM.png')

";
close($outhandle);
}
