#!/usr/bin/perl

use strict;
use warnings;
use Statistics::Regression;
use POSIX;
use Math::Amoeba qw( MinimiseND );
use Getopt::Long;
use Math::Matrix;
use Statistics::Basic qw(:all);

my $output_pair_separations=0;
my $data={};
my @labels=();
my $usage=join("\n",
	       "usage::Align_L1L2_L3.perl",
	       "Align_L1L2_L3.perl -mod=[nom|AB|CD|CD_220126|v3.12_RC]",
	       "\t[ --compare_sample=<sample_file> ][ --output_sample=<sample_file> ][ --gen_sens=<sensitivity_matrix>]",
	       "\t[--error=<measurement_error>] [--L2_only]");

my $deg=atan2(1,1)/45.0;
my $mod="nom";
my $sample_input;
my $sample_output;
my $gen_sensitivity_matrix;
my $meas_err=0.0;
my $only_L2=0;

GetOptions("mod=s" => \$mod,
	   "compare_sample=s" => \$sample_input,
	   "output_sample=s"  => \$sample_output,
	   "error=s"          => \$meas_err,
	   "L2_only"          => \$only_L2,
	   "gen_sens=s" => \$gen_sensitivity_matrix) || die "Error in commandline arguments! exiting..\n".$usage;

$mod=uc $mod;
$mod="nom" if ($mod eq "NOM");
$mod="ab"  if ($mod eq "AB");
$mod="cd"  if ($mod eq "CD");
$mod="cd_220126"  if ($mod eq "CD_220126");
$mod="v3.12_RC"  if ($mod eq "V3.12_RC");

die "error in commandline arguments! exiting..\n".$usage if (($mod ne "nom") &&
							     ($mod ne "ab") &&
							     ($mod ne "cd") &&
							     ($mod ne "cd_220126") &&
							     ($mod ne "v3.12_RC"));
my $perturb_L1S1={("nom"=>{"dir"=>1},
		   "ab" =>{("dir"=>1,
			    "tx" => -0.2875,"ty" => -0.2806,"tz" => +100e-3,
			    "rx" => -3.42e-3*$deg,"ry" => +6.71e-3*$deg)},
		   "cd" =>{("dir"=>1,
			    "tx" => -0.1204,"ty" => +0.1278,"tz" => +113.8e-3,
			    "rx" => -6.55e-3*$deg,"ry" => -1.15e-3*$deg)},
		   "cd_220126" => {("dir" => 1, # following values read off of ZEMAX optimization results, reversed signs of X & Z components: LSST_Ver_3.C_Baseline_Design_MOD_A_L1L2_asbuilt_L1L2hpod_FPmatch__M2hpod_fixed, referenced against the "nominal" axisymmetric, corrected version of v3.11: (LSST_Ver_3.11_Baseline_Design_L1S2-L2S1=412.571_Camera_AS):
				    "tx"  =>-0.137,
				    "ty"  =>-0.706,
				    "tz"  =>+0.120, # this is relative to the fiducial position 
				    "rx"  =>-0.011*$deg,
				    "ry"  =>-3.56e-3*$deg)},
		   "v3.12_RC_test" => {("dir" => 1,
				   "tz"  => +0.091, # relative to fiducial position
				   "tx"  => +0.037,
				   "ty"  => +0.120,
				   "rx"  => -6.115e-3*$deg,
				   "ry"  => -1.892e-3*$deg)},
		   "v3.12_RC" => {("dir" => 1,
				   "tz"  => +0.000, # say fiducial position is axisymmetric version in v3.12_RC
				   "tx"  => +0.007,
				   "ty"  => -0.022,
				   "rx"  => -6.68e-3*$deg,
				   "ry"  => -2.07e-3*$deg)},
    )};

printf STDERR "READING IN L3 SMR locations data..\n(12 L3 axial mounted SMRs wrt L3 OCF/DATUM A\n plus 12 L3 radial mounted SMRs wrt L3 OCF/DATUM A) .. ";

open(F,"<","L3_SMR_locations.dat") || die;

my ($tooling_ball_size_correction,$tooling_ball_diameter,$smr_nest_cone_incl_angle)=
    (1,12.000,30*$deg);

while (<F>) {
    next if (!/SPHERE/);
    chomp;
    my ($label,$x,$y,$z)=split(' ',$_);
    $data->{$label}={"X"=>$x,
		     "Y"=>$y,
		     "Z"=>$z};

    if (($label =~ /MEAS/) && $tooling_ball_size_correction) {
	my $correction=0.5*(12.700-$tooling_ball_diameter)/cos($smr_nest_cone_incl_angle);
	my $smr_id=(split("_",$label))[1];
	if ($smr_id<=12) {
	    # axial SMR, smr_ix ranges from 1 to 12:
	    $data->{$label}->{"Z"} += -1*$correction;
	} else {
	    # radial SMR, smr_id ranges from 13 to 24.
	    my $scalar = 1 + $correction/sqrt(pow($data->{$label}->{"X"},2)+
					      pow($data->{$label}->{"Y"},2));
	    $data->{$label}->{"X"} *= $scalar;
	    $data->{$label}->{"Y"} *= $scalar;
	}
    }
    push(@labels,$label);
}
close(F);
printf STDERR "done!\n";

printf STDERR "READING IN L1-L2 SMR LOCATIONS DATA:\n(6 L1 SMRs and 7 L2 SMRs in assembly\n plus 12 L1 cell measurements wrt L1S1 OCF\n plus 12 L2 Cell measurements wrt L2S2 OCF).. ";
open(F,"<","L1L2_SMR_locations.dat") || die;
while (<F>) {
    next if (!/SPHERE/);
    chomp;
    my ($label,$x,$y,$z)=split(' ',$_);
    $data->{$label}={"X"=>$x,
		     "Y"=>$y,
		     "Z"=>$z};
    push(@labels,$label);
}
close(F);
printf STDERR "done!\n\n";

printf STDERR "READING IN \"raw\" L1-L2 SMR LOCATIONS DATA:\n(6 L1 SMRs and 7 L2 SMRs in assembly\n plus 12 L1 cell measurements wrt L1S1 OCF\n plus 12 L2 Cell measurements wrt L2S2 OCF).. ";
open(F,"<","L1L2_SMR_NOM.dat") || die;
while (<F>) {
    next if (!/SPHERE/);
    chomp;
    my ($label,$x,$y,$z)=split(' ',$_);
    $data->{$label}={"X"=>$x,
		     "Y"=>$y,
		     "Z"=>$z};
    push(@labels,$label);
}
close(F);
printf STDERR "done!\n\n";

# now generate the point cloud for plotting

my @L3labels   = map { (($_ =~ /REF/) && ($_ =~ /THEO/))? ($_) : ()} @labels;
my @L3labels_meas   = map { (($_ =~ /REF/) && ($_ =~ /MEAS/))? ($_) : ()} @labels;
my @L1L2labels_nom = map { (($_ =~ /SMR/) && ($_ =~ /NOM/))? ($_) : ()} @labels;
my @L1L2labels = map { (($_ =~ /SMR/) && ($_ =~ /MEAS/))? ($_) : ()} @labels;

my @L1labels   = map { ($_ =~ /2560710/)? ($_) : ()} @labels;
my @L2labels   = map { ($_ =~ /2535589/)? ($_) : ()} @labels;

my @L1L2labels_oasm1 = map { (($_ =~ /SMR/) && ($_ =~ /OASM1/))? ($_) : ()} @labels;
my @L1L2labels_oasm2 = map { (($_ =~ /SMR/) && ($_ =~ /OASM2/))? ($_) : ()} @labels;
my @L1L2labels_oasm3 = map { (($_ =~ /SMR/) && ($_ =~ /OASM3/))? ($_) : ()} @labels;

printf STDERR "Filing away plottable representations of the datasets just read .. ";

save_measurements_file($data,[@L1L2labels],     "L1L2SMR_measurements.qdp");
save_measurements_file($data,[@L1L2labels_nom], "L1L2SMR_nominal.qdp", 1);
save_measurements_file($data,[@L1labels],       "L1SMR_BS2560710.qdp");
save_measurements_file($data,[@L2labels],       "L2SMR_BS2535589.qdp");
save_measurements_file($data,[@L3labels],       "L3SMR_measurements.qdp", 1);
save_measurements_file($data,[@L3labels_meas],  "L3SMR_measurements_meas.qdp", 1);

save_measurements_file($data,[@L1L2labels_oasm1],  "L1L2SMR_2565590_1.qdp", 1);
save_measurements_file($data,[@L1L2labels_oasm2],  "L1L2SMR_2565590_2.qdp", 1);
save_measurements_file($data,[@L1L2labels_oasm3],  "L1L2SMR_2565590_3.qdp", 1);

printf STDERR "done!\n\n";

if (1) {
    my $forward={"dir" => +1};
    # explain the table on 2565590 p.12
    my $tm_L1_p12  = match_pairs([@L1L2labels_oasm1],[@L1labels],["SMR\\d+","L1_Cell"]);
    my $dtf_L1_p12 = get_transform($data,$tm_L1_p12,$forward);
    my $tm_L2_p12  = match_pairs([@L1L2labels_oasm1],[@L2labels],["SMR\\d+","L2_Cell"]);
    my $dtf_L2_p12 = get_transform($data,$tm_L2_p12,$forward);
    # the following transformation query seems to answer the relationship between L1L2labels_oasm1 & L1L2labels..
    my $tm_L1L2_p12  = match_pairs([@L1L2labels_oasm1],[@L1L2labels],["SMR\\d+","L\\d_Cell"]);
    my $dtf_L1L2_p12;
    my $dtf_alcondL1_p19;
    my $dtf_alcondL2_p19;
    my $thingo;
    my $tm_L1_p15  = match_pairs([@L1L2labels_oasm1],[@L1labels],["SMR\\d+","L1_Cell"]);
    my $tm_L2_p15  = match_pairs([@L1L2labels_oasm1],[@L2labels],["SMR\\d+","L2_Cell"]);
    my $tm_L1_p19  = match_pairs([@L1L2labels_oasm3],[@L1labels],["SMR\\d+","L1_Cell"]);
    my $tm_L2_p19  = match_pairs([@L1L2labels_oasm3],[@L2labels],["SMR\\d+","L2_Cell"]);

    if (1) {
	$thingo={(%{$forward})};
	printf STDERR "%s\n13 SMRs of L1L2 (150-01520_p16) transformed to match corresponding positions of (2565590_p15)\n",join('',("#")x80);
	$dtf_L1L2_p12 = get_transform($data,$tm_L1L2_p12,$thingo);
	$thingo={(%{$forward},
		  "scale" => {("dir" => 0,	   "rx" => 0,	   "ry" => 0,
			       "rz" => 0,	   "tx" => 1,	   "ty" => 1,	   "tz" => 1)})};
	printf STDERR "%s\n13 SMRs of L1L2 (150-01520_p16) transformed to match corresponding positions of (2565590_p15) - trivial rotation\n",join('',("#")x80);
	$dtf_L1L2_p12 = get_transform($data,$tm_L1L2_p12,$thingo);
	
	$thingo={(%{$forward})};
	printf STDERR "%s\n6 SMRs of L1 (SER2560710) transformed to match corresponding positions of (2565590_p15)\n",join('',("#")x80);
	$dtf_L1L2_p12 = get_transform($data,$tm_L1_p15,$thingo);
	printf STDERR "%s\n7 SMRs of L2 (SER2535589) transformed to match corresponding positions of (2565590_p15)\n",join('',("#")x80);
	$dtf_L1L2_p12 = get_transform($data,$tm_L2_p15,$thingo);
	printf STDERR "%s\n",join('',("#")x80);

	$thingo={(%{$forward})};
	printf STDERR "%s\n12 SMRs of L1 (SER2560710) transformed to match corresponding positions of (2565590_p19)\n",join('',("#")x80);
	my $dtf_L1_MSRp19 = get_transform($data,$tm_L1_p19,$thingo);
	printf STDERR "%s\n12 SMRs of L2 (SER2535589) transformed to match corresponding positions of (2565590_p19)\n",join('',("#")x80);
	my $dtf_L2_MSRp19 = get_transform($data,$tm_L2_p19,$thingo);
	printf STDERR "%s\n",join('',("#")x80);

	# ********************** OLD DIAGNOSIS ***************************
	# we suspect there was a sign error committed when the authors generated MSRp19. Here we can produce a transformation that puts things back, and find the transformation needed to fit the lens-level SMR positions to the (corrected) configuration.

	#	my $transform_to_correct={"dir" => 1,
	#				      "tx" => 2*(-0.061), "ty" => 2*(+0.03), "tz" => -464.069}; # cf. -958.95 mm ± 0.2 mm relative to the interface datum on its support structure defined in the L1-L2 Assembly Interface Definition Drawing (LCA-77; Note 3 not applicable) (OPT-L1L2STRUC019).
	#	# ** Difference between 464.069 and 958.950 is 494.881, which is the separation between L1S1 & L2S2. **
	#	*************** NB - not reporting SMR positions relative to the Interface DATUM ******************
	#                       ________________________________________________________________	
	
	# in the L1L2 assembly optical coordinate system if the sign error was committed as suspected. Origin is the intersection of L2S1 plane with the aspheric axis of L2S2.
	# my $transform_to_correct={"dir" => 1, "tx" => 2*(-0.061), "ty" => 2*(+0.03), "tz" => 0};

	# ********************** CURRENT DIAGNOSIS ***************************
	my $transform_to_correct={"dir" => 1, "tx" => 1*(-0.061), "ty" => 1*(+0.03), "tz" => 0};  # THIS WOULD BRING L2S2 APEX ALIGNED WITH WITH OPTICAL FRAME AXIS
	
	# ********************** HYPER-CURRENT DIAGNOSIS ***************************
	# the asphere decenter is not coincident with the thinnest part of L2 optic. With wedge consistent with zero, measuring relative to the
	# L2S2 Asphere Axis was the wrong thing to do.. Koby explains that the misaligned asphere imparts a small coma but that's not measurable in CGH test (only tilts are)
	# so OASM3 is considered to have the origin coincident with the CGH beam (with overall tilt removed). SO WE CAN FORGET ABOUT OASM4 AND USE OASM3 only.

	if (0) {
	    my @L1L2labels_oasm4 = @L1L2labels_oasm3;
	    # and alter their labels.
	    map { $_ =~ s/OASM3/OASM4/ } @L1L2labels_oasm4;
	    # populate the new labels:
	    populate_transformed_coords($data,[@L1L2labels_oasm4],
					$transform_to_correct,
					[@L1L2labels_oasm3]);
	    # prepare the pairs:
	    my $tm_L1_p19_corrected  = match_pairs([@L1L2labels_oasm4],[@L1labels],["SMR\\d+","L1_Cell"]);
	    my $tm_L2_p19_corrected  = match_pairs([@L1L2labels_oasm4],[@L2labels],["SMR\\d+","L2_Cell"]);
	    
	    # transformation done, now need to fit:
	    printf STDERR "%s\n",join('',("#")x80);
	    printf STDERR "%s\n12 SMRs of L1 (SER2560710) transformed to match *corrected* positions of (2565590_p19)\n",join('',("#")x80);
	    my $dtf_L1_MSRp19_corrected = get_transform($data,$tm_L1_p19_corrected,$forward);
	    printf STDERR "%s\n12 SMRs of L2 (SER2535589) transformed to match *corrected* positions of (2565590_p19)\n",join('',("#")x80);
	    my $dtf_L2_MSRp19_corrected = get_transform($data,$tm_L2_p19_corrected,$forward);
	    printf STDERR "%s\n",join('',("#")x80);
	    
	    # save this pattern of SMRs enven though it's derived from the tabulated numbers of p.19..
	    save_measurements_file($data,[@L1L2labels_oasm4],  "L1L2SMR_2565590_4_p19corrected.qdp", 1);
	    
	    # The L1L2labels_oasm4 set is now centered with respect to the optical axis
	    # and referenced against the L2S1 apex. If it were not for the non-axisymmetric details, all we would need to
	    # do at this point would be to translate the L3labels_meas set along z and combine with OASM4.

	    # now for L3.
	}
	# REMOVE BEFORE FLIGHT   exit;
    }
    # the L1-L2 part is complete, @L1L2labels_oasm3 are the labels, with coordinates reported relative to the Assembly Optical Coordinate System which coincides with the CGH beam axis under NULL condition (or close to it)
    # in oasm3 set, the origin is the "apex" of L2S1.
    # move on to handle L3..
    {
	# do this (L3 work) without using triads. Too confusing.
	# the set of L3 labels relative to L3IF is called @L3labels_meas
	my @L3labels_oasm3=@L3labels_meas;
	map {$_ =~ s/^/L3_/;$_ =~ s/A3/OASM3/} @L3labels_oasm3;
	# printf STDERR "%s\n",join("\n","L3LABELS_OASM3:",@L3labels_oasm3);
	# transform to show these together.
	my $L3_tf = {("dir" => -1,
		      "tx" => 0.0864, # Position of L3S1 apex wrt frame "MECA" system, transformed by rotation to CCS system (where SMR coords are reported)
		      "ty" => 0.1667, # Position of L3S1 apex wrt frame "MECA" system, transformed by rotation to CCS system (where SMR coords are reported)
		      "tz" => (+3.27 # L3S2 Apex minus L3IF
			       + sum(-30.05,-346.237,0,-17.90,-54.1,-60)) # SUM of thicknesses of L2, L2 to F, F, F to L3, L3. (CHECK THESE NUMBERS -- see fix below.)
	    )};
	# Kludge for now, using LCN-1898 numbers:
	$L3_tf->{"tz"} = -505.31; # different from the (initial) ZEMAX stackup by *** 293 microns ***

############################################################################################################
	# "nominal" axismmetric position based on a corrected version of v3.11 (LSST_Ver_3.11_Baseline_Design_L1S2-L2S1=412.571_Camera_AS): 
	$L3_tf->{"tz"} = (+3.27 # L3S2 Apex minus L3IF
			  + sum(-30.05,-346.248,0,-17.90,-54.1,-60)); # SUM of thicknesses of L2, L2 to F, F, F to L3, L3. (using LSST_Ver_3.11_Baseline_Design_L1S2-L2S1=412.571_Camera_AS):
	# L3_tf->{"tz"} for the fiducial (axisymmetric) position works out to -505.028mm, different from the v3.11 baseline by 11um and different from the LCN-1898 "implied baseline" by 282um.

############################################################################################################

############################################################################################################
	# axial position based on v3.12_RC (LSST_Ver_3.12_Baseline_Design_RC):
	if ($mod eq "v3.12_RC") {
	    printf "*************************** ALIGNMENT MODEL: v3.12_RC ******************************\n";
	    $L3_tf->{"tz"} = (+3.27 # L3S2 Apex minus L3IF
			      + sum(-30.05,-346.237,0,-17.90,-54.1,-60)); # SUM of thicknesses of L2, L2 to F, F, F to L3, L3. (using LSST_Ver_3.12_RC):
	}
	# L3_tf->{"tz"} for the fiducial (axisymmetric) position works out to -505.017mm, different from the LCN-1898 value (-505.31mm) by 293um.

############################################################################################################
	
	populate_transformed_coords($data,[@L3labels_oasm3],[$L3_tf],[@L3labels_meas]);
	save_measurements_file($data,[@L1L2labels_oasm3,@L3labels_oasm3],"L1L2L3_oasm3.qdp",1);
	if (1) {
	    open(SAMPLE,">","L3_SMRs.txt") || die;
#	    printf SAMPLE "# UNCOMMENT AND EDIT THE FOLLOWING ACCORDING TO INDIVIDUAL L1L2 SMR MEASUREMENTS\n";
#	    printf SAMPLE "# AFTER L3 SMR POSITIONS HAVE BEEN USED TO DEFINE THE COORDINATE SYSTEM:\n";
	    printf SAMPLE "# LABEL X Y Z\n";
	    foreach my $label (sort {(split("_",$a))[2] <=> (split("_",$b))[2]} @L3labels_oasm3) {
		printf SAMPLE "%s%-40s %10s %10s %10s\n",($label =~ /Cell/)?("# "):(""),$label,@{$data->{$label}}{"X","Y","Z"};
	    }
	}
	# but wait, there's more.
	# The nominal ends of the struts are captured in LCA-77. They are called STRUT_A{1,2,3,4,5,6}. Axial separation between planes is 696.831mm.
	# -Z end has STRUT_A[12] joining at a point, STRUT_A[34] joining, STRUT_A[56] joining, etc. The intersection points are equally spaced (120 deg) on diameter of 1504.00mm.
	# +Z end has STRUT_A[12] and STRUT_A[34] each separated by 34.39 deg (struts closer to (+/-)X-axis are 9.07 deg from the (+/-)X-axis),
	# and STRUT_A[56] separated by 48.655 deg - all on a bolt circle with diameter of 1590.00mm.
	my $r_isect=1504.00/2.0;
	my $r_foot_bc=1590.00/2.0;
	my $theta_foot={("A1" => 9.07*$deg,	 "A2" => (9.07+34.39)*$deg,	"A3" => (180-(9.07+34.39))*$deg,
			 "A4" => (180-9.07)*$deg,"A5" => (-90-48.655/2.0)*$deg,"A6" => (-90+48.655/2.0)*$deg)};
	my $theta_isect={("A1" => (-90+120)*$deg,"A2" => (-90+120)*$deg,	"A3" => (-90+2*120)*$deg,
			  "A4" => (-90+2*120)*$deg,"A5" => (-90+3*120)*$deg,"A6" => (-90+3*120)*$deg )};
	my $l1s1_l2s1_delta_z = -494.881; # this is measured relative to L2S1: L1S1 Vertex is -494.881mm in z from L2S1
	my $c_if_l1s1 = +959.321; # L1S1 is -959.321mm in z from Camera Interface Plane.
	my $camera_interface_z = $c_if_l1s1 + $l1s1_l2s1_delta_z; # straight sum of the above
	
	my $struts={("A1_FOOT"  => {("X" => $r_foot_bc*cos($theta_foot->{"A1"}),   "Y" => $r_foot_bc*sin($theta_foot->{"A1"}),    "Z" => $camera_interface_z)},
		     "A2_FOOT"  => {("X" => $r_foot_bc*cos($theta_foot->{"A2"}),   "Y" => $r_foot_bc*sin($theta_foot->{"A2"}),    "Z" => $camera_interface_z)},
		     "A3_FOOT"  => {("X" => $r_foot_bc*cos($theta_foot->{"A3"}),   "Y" => $r_foot_bc*sin($theta_foot->{"A3"}),    "Z" => $camera_interface_z)},
		     "A4_FOOT"  => {("X" => $r_foot_bc*cos($theta_foot->{"A4"}),   "Y" => $r_foot_bc*sin($theta_foot->{"A4"}),    "Z" => $camera_interface_z)},
		     "A5_FOOT"  => {("X" => $r_foot_bc*cos($theta_foot->{"A5"}),   "Y" => $r_foot_bc*sin($theta_foot->{"A5"}),    "Z" => $camera_interface_z)},
		     "A6_FOOT"  => {("X" => $r_foot_bc*cos($theta_foot->{"A6"}),   "Y" => $r_foot_bc*sin($theta_foot->{"A6"}),    "Z" => $camera_interface_z)},
		     "A1_ISECT" => {("X" => $r_isect*cos($theta_isect->{"A1"}),    "Y" => $r_isect*sin($theta_isect->{"A1"}),     "Z" => $camera_interface_z-696.831)},
		     "A2_ISECT" => {("X" => $r_isect*cos($theta_isect->{"A2"}),    "Y" => $r_isect*sin($theta_isect->{"A2"}),     "Z" => $camera_interface_z-696.831)},
		     "A3_ISECT" => {("X" => $r_isect*cos($theta_isect->{"A3"}),    "Y" => $r_isect*sin($theta_isect->{"A3"}),     "Z" => $camera_interface_z-696.831)},
		     "A4_ISECT" => {("X" => $r_isect*cos($theta_isect->{"A4"}),    "Y" => $r_isect*sin($theta_isect->{"A4"}),     "Z" => $camera_interface_z-696.831)},
		     "A5_ISECT" => {("X" => $r_isect*cos($theta_isect->{"A5"}),    "Y" => $r_isect*sin($theta_isect->{"A5"}),     "Z" => $camera_interface_z-696.831)},
		     "A6_ISECT" => {("X" => $r_isect*cos($theta_isect->{"A6"}),    "Y" => $r_isect*sin($theta_isect->{"A6"}),	  "Z" => $camera_interface_z-696.831)})};
	my @strut_names=();
	foreach my $point (sort glob("A{1,2,3,4,5,6}_{FOOT,ISECT}")) {
	    push(@strut_names,$point);
	    $data->{$point}=$struts->{$point};
	}
	# now push these labels onto the lists: A[1-6]_FOOT onto @L3labels_oasm3 and A[1-6]_ISECT onto @L1L2labels_oasm3.
	save_measurements_file($data,[@L1L2labels_oasm3,@L3labels_oasm3],"L1L2L3_oasm3_nostruts.qdp",0);

	if (defined($sample_input)) { # read in sample.txt
	    # read in "sample.txt" L2 SMR positions to see how they jive with the expected oasm3 positions
	    my $input=$sample_input;
	    my @sample_labels=();
	    printf STDERR "reading data from file: $input\n";
	    open(F,"<",$input) || die;
	    while (<F>) {
		next if (/^#/);
		next if (!/SPHERE/);
		chomp;
		my ($label,$x,$y,$z)=split(' ',$_);
		$data->{$label}={"X"=>$x,"Y"=>$y,"Z"=>$z};

		# won't be necessary from now on, this was only for first adjustment influenced by incorrect targets

		my $one_time_kludge=0; 
		if ($one_time_kludge && $tooling_ball_size_correction) {
		    my $correction=0.5*(12.700-$tooling_ball_diameter)/cos($smr_nest_cone_incl_angle);
		    $data->{$label}->{"Z"} -= $correction;
		}

		push(@sample_labels,$label);
	    }
	    close(F);
	    # prepare @L1L2labels_oasm3_rota (perturbed target positions)
	    my @L1L2labels_oasm3_rota=@L1L2labels_oasm3;
	    map {$_ =~ s/$/_rota/} @L1L2labels_oasm3_rota;
	    populate_transformed_coords($data,
					[@L1L2labels_oasm3_rota],
					[({("dir" => -1,"tz" => $l1s1_l2s1_delta_z)}, # rotation about L1S1 vertex (nom)
					  $perturb_L1S1->{$mod},
					  {("dir" => +1,"tz"=>  $l1s1_l2s1_delta_z)})],
					[@L1L2labels_oasm3]);
	    
	    # now compare the SAMPLE labels with @L1L2labels_oasm3_rota:
	    my $sample_pairs=[];
	    # get the pair ordering right:
	    push(@{$sample_pairs},@{match_pairs([@sample_labels],[@L1L2labels_oasm3_rota],["L2_SMR\\d+"])});
	    push(@{$sample_pairs},@{match_pairs([@sample_labels],[@L1L2labels_oasm3_rota],["L2_Cell_SMR\\d+"])});
	    
	    # move on to fit.
	    my $tm_init={%{$forward}};
	    my $fit_tf=get_transform($data,$sample_pairs,$tm_init);
	    # now include the strut coordinates into the dataset and figure out the adjustment.
	    # (have @sample_labels, want @L1L2labels_oasm3_rota):
	    my @target_foot=glob("A{1,2,3,4,5,6}_FOOT");
	    my @target_isect=glob("A{1,2,3,4,5,6}_ISECT");
	    my @target_isect_sample=glob("A{1,2,3,4,5,6}_ISECT_SAMPLE");
	    populate_transformed_coords($data,[@target_isect_sample],[$fit_tf],[@target_isect]);
	    my $strutlengths_nominal = {
		map { "A".$_ =>
			  scalar_separation($data->{"A".$_."_FOOT"},
					    $data->{"A".$_."_ISECT"}) } ( 1..6 ) };
	    my $strutlengths_sample = {
		map { "A".$_ => scalar_separation($data->{"A".$_."_FOOT"},
						  $data->{"A".$_."_ISECT_SAMPLE"}) } ( 1..6 ) };

	    printf STDERR "\nL1-L2 Assy rel. L3.. %s:\n",join(" ",map {sprintf("%s => %s",$_,$fit_tf->{$_})} sort keys %{$fit_tf});
	    map {printf STDERR "A".$_."delta_length = %f\n",
		     -1*($strutlengths_sample->{"A".$_}-$strutlengths_nominal->{"A".$_}) } (1..6);
	    exit;
	}
	
	my @model_tf_oasm3;
	
	if (defined($sample_output)) { # output sample.txt
	    # this snippet outputs a sample datafile to edit with actual L1L2 SMR coordinates reported relative to L3 SMRs (not checking for L3 pattern fit):
	    # first provide a little bit of an inadvertent transformation:
	    my $err_tf = {("dir"=>1,
			   #			   "rx"=>50e-6,
			   "tx"=>1)};
	    my @L1L2labels_oasm3_tf=@L1L2labels_oasm3;
	    map {$_ =~ s/$/_ACTUAL/} @L1L2labels_oasm3_tf;

	    populate_transformed_coords($data,
					[@L1L2labels_oasm3_tf],
					[({("dir" => -1,"tz" => $l1s1_l2s1_delta_z)}, # rotation about L1S1 vertex (nom)
					  $perturb_L1S1->{$mod},
					  {("dir" => +1,"tz"=>  $l1s1_l2s1_delta_z)})],
					[@L1L2labels_oasm3]);

	    @model_tf_oasm3=@L1L2labels_oasm3_tf;
	    
	    # add in some measurement error
	    foreach my $pt (@L1L2labels_oasm3_tf) {
		foreach my $ax ("X","Y","Z") {
		    $data->{$pt}->{$ax} += $meas_err*(2*rand()-1);
		}
	    }
	    
	    open(SAMPLE,">",$sample_output) || die;
	    printf SAMPLE "# UNCOMMENT AND EDIT THE FOLLOWING ACCORDING TO INDIVIDUAL L1L2 SMR MEASUREMENTS\n";
	    printf SAMPLE "# AFTER L3 SMR POSITIONS HAVE BEEN USED TO DEFINE THE COORDINATE SYSTEM:\n";
	    printf SAMPLE "# LABEL X Y Z\n";
	    if ($only_L2) {
		foreach my $label (sort {(split("_",$a))[1] <=> (split("_",$b))[1]} @L1L2labels_oasm3_tf) {
		    next if ($label !~ /L2/);
		    next if ($label =~ /Cell/);
		    printf SAMPLE "%s%-40s %10s %10s %10s\n",($label =~ /Cell/)?("# "):(""),$label,@{$data->{$label}}{"X","Y","Z"};
		}
	    } else {
		foreach my $label (sort {(split("_",$a))[1] <=> (split("_",$b))[1]} @L1L2labels_oasm3_tf) {
		    next if ($label !~ /AR/);
		    printf SAMPLE "%s%-40s %10s %10s %10s\n","",$label,@{$data->{$label}}{"X","Y","Z"};
		}
		foreach my $label (sort {(split("_",$a))[1] <=> (split("_",$b))[1]} @L1L2labels_oasm3_tf) {
		    next if ($label !~ /L1/);
		    printf SAMPLE "%s%-40s %10s %10s %10s\n","",$label,@{$data->{$label}}{"X","Y","Z"};
		}
		foreach my $label (sort {(split("_",$a))[1] <=> (split("_",$b))[1]} @L1L2labels_oasm3_tf) {
		    next if ($label !~ /L2/);
		    next if ($label =~ /Cell/);
		    printf SAMPLE "%s%-40s %10s %10s %10s\n","",$label,@{$data->{$label}}{"X","Y","Z"};
		}
		foreach my $label (sort {(split("_",$a))[1] <=> (split("_",$b))[1]} @L1L2labels_oasm3_tf) {
		    next if ($label !~ /L2/);
		    next if ($label !~ /Cell/);
		    printf SAMPLE "%s%-40s %10s %10s %10s\n","",$label,@{$data->{$label}}{"X","Y","Z"};
		}
	    }
	    close(SAMPLE);
	}

	push(@L1L2labels_oasm3,glob("A{1,2,3,4,5,6}_ISECT"));
	push(@L3labels_oasm3,glob("A{1,2,3,4,5,6}_FOOT"));
	save_measurements_file($data,[@L1L2labels_oasm3,@L3labels_oasm3],"L1L2L3_oasm3_struts.qdp",0);

	# record strut lengths
	my $strutlengths_nominal = {
	    map { "A".$_ =>
		      scalar_separation($data->{"A".$_."_FOOT"},
					$data->{"A".$_."_ISECT"}) } ( 1..6 ) };

	# test 1. copy & transform the @L1L2labels_oasm3 (which now include struts)
	my @L1L2labels_oasm3_rota = @L1L2labels_oasm3;

	map {$_ =~ s/$/_rota/} @L1L2labels_oasm3_rota;
	if (defined($gen_sensitivity_matrix)) { # generate sensitivity matrix
	    my $SENSMAT = [];
	    open(SENSMAT,">",$gen_sensitivity_matrix) || die;
	    foreach my $single_elem_tf ({"dir"=>+1,"tx"=>1.0},
					{"dir"=>+1,"ty"=>1.0},
					{"dir"=>+1,"tz"=>1.0},
					{"dir"=>+1,"rx"=>100e-6},
					{"dir"=>+1,"ry"=>100e-6},
					{"dir"=>+1,"rz"=>100e-6}) {
		
		populate_transformed_coords($data,
					    [@L1L2labels_oasm3_rota],
					    [({("dir" => -1,"tz" => $l1s1_l2s1_delta_z)}, # rotation about L1S1 vertex (nom)
					      $single_elem_tf,
					      {("dir" => +1,"tz"=>  $l1s1_l2s1_delta_z)})],
					    [ @L1L2labels_oasm3 ]); # would like to use @L1L2labels_oasm3_tf instead of @L1L2labels_oasm3 (to accommodate different models)
		
		my $strutlengths_rota = { map { "A".$_ => scalar_separation($data->{"A".$_."_FOOT"},$data->{"A".$_."_ISECT_rota"}) } ( 1..6 ) };

		printf SENSMAT "\nL1S1Vertex rel. L3.. %s:\n",join(" ",map {sprintf("%s => %s",$_,$single_elem_tf->{$_})} sort keys %{$single_elem_tf});
		map {printf SENSMAT "A".$_."delta_length = %f\n",
			 $strutlengths_rota->{"A".$_}-$strutlengths_nominal->{"A".$_}} (1..6);
		#	    save_measurements_file($data,[@L1L2labels_oasm3_rota,@L3labels_oasm3],"L1L2L3_oasm3_struts_rota.qdp",0);
		# pack into matrix.
		push (@{$SENSMAT},[map {$strutlengths_rota->{"A".$_}-$strutlengths_nominal->{"A".$_}} (1..6)]);
	    }
	    close(SENSMAT);
	    # now invert/compute the sensitivity matrix etc.
	    my $A=Math::Matrix->new(@{$SENSMAT});
	    $A=$A->transpose();
	    $A->print("A\n");
	    my $B=$A->inv;
	    $B->print("B\n");
	    my $trials=1000000;
	    my $results = [];
	    for (my $i=0;$i<$trials;$i++) {
		my $errs=[[]];
		my $amp=0.1;
		map {$errs->[0]->[$_]=$amp*(0.5-rand(1.0))} (0..5);
		my $C=Math::Matrix->new_from_cols($errs);
		my $D=$B->mul($C);
#		$D->print("D\n");
		$D=$D->transpose;
		$results->[$i]=[@{$D->[0]}];
#		$D->print("D_transpose\n");
#		printf "%s\n",join(" ",map {sprintf("%f",$_)} @{$D->[0]});
	    }
	    foreach my $entry (0..2) {
		my @list=map {($results->[$_]->[$entry])} (0..$trials-1);
		printf "STATS FOR ENTRY $entry: mean %10f stddev %10f\n",mean(@list),stddev(@list);
	    }
	    foreach my $entry (3..5) {
		my @list=map {($results->[$_]->[$entry])} (0..$trials-1);
		printf "STATS FOR ENTRY $entry(deg): mean %10f stddev %10f\n",57.3e-4*mean(@list),57.3e-4*stddev(@list);
	    }
	}

	if (1) { # generate offax_tuning settings
	    open(NAS,">","nonaxisymmetric_settings.txt") || die;
	    my $nas_list={("NOMINAL" => {"dir" => 1},
			   "Solution_A_B rel. NOM" => $perturb_L1S1->{"ab"},
			   "Solution_C_D rel. NOM" => $perturb_L1S1->{"cd"},
			   "Solution_C_D_220126 rel. NOM" => $perturb_L1S1->{"cd_220126"},
		)};
	    if ($mod eq "v3.12_RC") {
		$nas_list={("NOMINAL" => {"dir" => 1},
			    "Solution_v3.12_RC rel. NOM" => $perturb_L1S1->{"v3.12_RC"})};
	    }
	    foreach my $nas_soln (sort keys %{$nas_list}) {
		my $single_elem_tf=$nas_list->{$nas_soln};
		populate_transformed_coords(
		    $data,
		    [@L1L2labels_oasm3_rota],
		    [({("dir" => -1,"tz" => $l1s1_l2s1_delta_z)}, # rotation about L1S1 vertex (nom)
		      $single_elem_tf,
		      {("dir" => +1,"tz"=>  $l1s1_l2s1_delta_z)})],
		    [@L1L2labels_oasm3]);
		
		my $strutlengths_rota = {
		    map { "A".$_ =>
			      scalar_separation($data->{"A".$_."_FOOT"},
						$data->{"A".$_."_ISECT_rota"}) } ( 1..6 ) };
		printf NAS "\nL1S1Vertex rel. L3 (%s).. %s:\n",$nas_soln,join(" ",map {sprintf("%s => %s",$_,$single_elem_tf->{$_})} sort keys %{$single_elem_tf});
		map {printf NAS "A".$_."delta_length = %f\n",
			 $strutlengths_rota->{"A".$_}-$strutlengths_nominal->{"A".$_}} (1..6);
		#	    save_measurements_file($data,[@L1L2labels_oasm3_rota,@L3labels_oasm3],"L1L2L3_oasm3_struts_rota.qdp",0);
	    }
	    close(NAS);
	}
	save_measurements_file($data,[@strut_names],"strut_points.qdp",0);
    }
    exit;
    
    $thingo={(%{$forward})};

    $dtf_L1L2_p12 = get_transform($data,$tm_L1L2_p12,$thingo);
    printf STDERR "%s\n",join('',("#")x80);

    $thingo={(%{$forward})};
    my $tm_alcond_oasm3_L1  = match_pairs([@L1L2labels_oasm3],[@L1labels],["SMR\\d+","L1_Cell"]);
    my $tm_alcond_oasm3_L2  = match_pairs([@L1L2labels_oasm3],[@L2labels],["SMR\\d+","L2_Cell"]);

    printf STDERR "%s\n() vs () (L1)\n",join('',("#")x80);
    $dtf_alcondL1_p19 = get_transform($data,$tm_alcond_oasm3_L1,$thingo);
    printf STDERR "%s\nOASM2 vs OASM3 (L2)\n",join('',("#")x80);
    $dtf_alcondL2_p19 = get_transform($data,$tm_alcond_oasm3_L2,$thingo);
    printf STDERR "%s\n",join('',("#")x80);

#    exit; COMMENTED OUT
    $thingo={(%{$forward},
	      "scale" => {("dir" => 0,	   "rx" => 0,	   "ry" => 0,
			   "rz" => 0,	   "tx" => 1,	   "ty" => 1,	   "tz" => 1)})};
    printf STDERR "%s\n13 SMRs of L1L2 (150-01520_p16) transformed to match corresponding positions of (2565590_p15) - trivial rotation\n",join('',("#")x80);
    $dtf_L1L2_p12 = get_transform($data,$tm_L1L2_p12,$thingo);
    printf STDERR "%s\n",join('',("#")x80);
    

#    exit; COMMENTED OUT
    # try looking for matches between L1L2labels_oasm1 and L1L2labels
    my $testmatch_L1_1_v_2 = match_pairs([@L1L2labels_oasm1],[@L1L2labels_oasm2],["SMR\\d+","L1_Cell"]);
    my $testmatch_L2_1_v_2 = match_pairs([@L1L2labels_oasm1],[@L1L2labels_oasm2],["SMR\\d+","L2_Cell"]);
    my $testmatch_L1_1_v_3 = match_pairs([@L1L2labels_oasm1],[@L1L2labels_oasm3],["SMR\\d+","L1_Cell"]);
    my $testmatch_L2_1_v_3 = match_pairs([@L1L2labels_oasm1],[@L1L2labels_oasm3],["SMR\\d+","L2_Cell"]);

    my $testmatch_L1_1 = match_pairs([@L1labels],[@L1L2labels_oasm1],["SMR\\d+","L\\d_Cell"]);
    my $testmatch_L1_2 = match_pairs([@L1labels],[@L1L2labels_oasm2],["SMR\\d+","L\\d_Cell"]);
    my $testmatch_L1_3 = match_pairs([@L1labels],[@L1L2labels_oasm3],["SMR\\d+","L\\d_Cell"]);
    my $testmatch_L2_1 = match_pairs([@L2labels],[@L1L2labels_oasm1],["SMR\\d+","L\\d_Cell"]);
    my $testmatch_L2_2 = match_pairs([@L2labels],[@L1L2labels_oasm2],["SMR\\d+","L\\d_Cell"]);
    my $testmatch_L2_3 = match_pairs([@L2labels],[@L1L2labels_oasm3],["SMR\\d+","L\\d_Cell"]);
    my $testmatch_3_L2 = match_pairs([@L1L2labels_oasm3],[@L2labels],["SMR\\d+","L\\d_Cell"]);
    #    printf "testmatch: %s\n",join("\n",sort map {join(",",@{$_})} @{$testmatch});

    printf STDERR "%s\nOASM1 to L1\n",join('',("#")x80);
    my $dtf_L1_1 = get_transform($data,$testmatch_L1_1,$forward);
    printf STDERR "%s\nOASM2 to L1\n",join('',("#")x80);
    my $dtf_L1_2 = get_transform($data,$testmatch_L1_2,$forward);
    printf STDERR "%s\nOASM3 to L1\n",join('',("#")x80);
    my $dtf_L1_3 = get_transform($data,$testmatch_L1_3,$forward);
    printf STDERR "%s\nOASM1 to L2\n",join('',("#")x80);
    my $dtf_L2_1 = get_transform($data,$testmatch_L2_1,$forward);
    printf STDERR "%s\nOASM2 to L2\n",join('',("#")x80);
    my $dtf_L2_2 = get_transform($data,$testmatch_L2_2,$forward);
    printf STDERR "%s\nOASM3 to L2\n",join('',("#")x80);
    my $dtf_L2_3 = get_transform($data,$testmatch_L2_3,$forward);
    printf STDERR "%s\nSER2535589 to 2565590_p13\n",join('',("#")x80);
    my $dtf_3_L2 = get_transform($data,$testmatch_3_L2,$forward);
    printf STDERR "%s\nSER2535589 to 2565590_p13 (trivial rotation)\n",join('',("#")x80);
    $thingo={(%{$forward},
	      "scale" => {("dir" => 0,	   "rx" => 0,	   "ry" => 0,
			   "rz" => 0,	   "tx" => 1,	   "ty" => 1,	   "tz" => 1)})};
    $dtf_3_L2 = get_transform($data,$testmatch_3_L2,$thingo);

    printf STDERR "%s\nL1(assy) OASM2 to OASM1\n",join('',("#")x80);
    my $dtf_L1_1_v_2 = get_transform($data,$testmatch_L1_1_v_2,$forward);
    printf STDERR "%s\nL2(assy) OASM2 to OASM1\n",join('',("#")x80);
    my $dtf_L2_1_v_2 = get_transform($data,$testmatch_L2_1_v_2,$forward);

    printf STDERR "%s\nL1(assy) OASM3 to OASM1\n",join('',("#")x80);
    my $dtf_L1_1_v_3 = get_transform($data,$testmatch_L1_1_v_3,$forward);
    printf STDERR "%s\nL2(assy) OASM3 to OASM1\n",join('',("#")x80);
    my $dtf_L2_1_v_3 = get_transform($data,$testmatch_L2_1_v_3,$forward);
    printf STDERR "%s\n",join('',("#")x80);
# exit; COMMENTED OUT
}
printf STDERR "IDENTIFYING SMR COORDINATES ACROSS DATA SETS .. \n";
# find transformation that joins the L1L2labels &  L1labels sets 

printf STDERR " .. matches for L1 labels:\n";
my $l1p = match_pairs([@L1L2labels],[@L1labels],["SMR\\d+","L\\d_Cell"]);
# my $l1l2ppp = match_pairs([@L1L2labels],[@L1L2labels_nom],["SMR\\d+","L\\d_Cell"]);
my $l1_in_l1l2ppp = match_pairs([@L1L2labels],[@L1L2labels_nom],["SMR\\d+","L1_Cell"]);

printf STDERR " .. matches for L2 labels:\n";
my $l2p = match_pairs([@L1L2labels],[@L2labels],["SMR\\d+","L\\d_Cell"]);
# the coordinate transformation between the above is a straight translation along Z of 524.922mm
my $l2_in_l1l2ppp = match_pairs([@L1L2labels],[@L1L2labels_nom],["SMR\\d+","L2_Cell"]);

printf STDERR "DONE!\n\n";
if (0) {
    my $null_tran={"rx" => 0,"ry" => 0,"rz" => 0,
		       "tx" => 0,"ty" => 0,"tz" => 0,"dir" => 1};
    printf STDERR "convar: %f\n",sqrt(constellation_var($data,$l1p,$null_tran));
    printf STDERR "convar: %f\n",sqrt(constellation_var($data,$l2p,$null_tran));
}

my $derived_transform={};

my $fixed_transform={("tx" => 0,		 "ty" => 0,		 "tz" => 0,
		      "rx" => 0.0,		 "ry" => 0.0,		 "rz" => 0.0,
		      "dir" => 0)};

my $translation_only={("tx" => 1,		 "ty" => 1,		 "tz" => 1,
		       "rx" => 0.0,		 "ry" => 0.0,		 "rz" => 0.0,
		       "dir" => 0)};

my $forward={"dir"=>+1};
my $forward_translation_only={("dir"=>+1,
			       "scale" => $translation_only)};
#$forward_translation_only=$forward;
$forward_translation_only={("dir"=>+1,
			    "scale" => $fixed_transform)};
my $backward={"dir"=>+1};
my $backward_translation_only={("dir"=>-1,
				"scale" => $translation_only)};

my $dt;

printf STDERR "\nFitting the nominal L1 SMRs to the actual (as measured) L1-L2 SMR configuration..\n";
#$dt=get_transform($data,$l1_in_l1l2ppp,$forward);
#printf STDERR "convar: %f\n",sqrt(constellation_var($data,$l1_in_l1l2ppp,$dt,"L1_smr_resids.qdp"));
$dt=get_transform($data,$l1p,$forward);
printf STDERR "convar: %f\n",sqrt(constellation_var($data,$l1p,$dt,"L1_smr_resids.qdp"));

printf STDERR "\nFitting the nominal L2 SMRs to the actual (as measured) L1-L2 SMR configuration..\n";
#$dt=get_transform($data,$l2_in_l1l2ppp,$forward);
#printf STDERR "convar: %f\n",sqrt(constellation_var($data,$l2_in_l1l2ppp,$dt,"L2_smr_resids.qdp"));
$dt=get_transform($data,$l2p,$forward);
printf STDERR "convar: %f\n",sqrt(constellation_var($data,$l2p,$dt,"L2_smr_resids.qdp"));


# fill out the full complement of L1 & L2 SMRs according to the transforms

printf STDERR "\nThis will be L1_rel_L1L2:\n";
$derived_transform->{"L1_rel_L1L2"}=get_transform($data,$l1p,$forward);

# copy over labels into a new list
my @L1_in_L1L2_labels=@L1labels;
# transform labels in place
map {($_ =~ s/BS2560710/L1L2assy_BS2560710/)?($_):()} @L1_in_L1L2_labels;
# populate the coordinates of the new list according to the transform
populate_transformed_coords($data,[@L1_in_L1L2_labels],
			    $derived_transform->{"L1_rel_L1L2"},[@L1labels]);

# print out the rms 3-d variance between the transformed SMR group wrt the original SMR group.
printf STDERR "check: N = %d; y = %f\n",scalar(@L1_in_L1L2_labels),sqrt(constellation_var($data,mux_into_pairs([@L1_in_L1L2_labels],[@L1labels])));


if (0) {
    # verify that the same transform can be extracted going forward with the group order as
    # can be extracted going backward with the reversed group order: (this works, by the way)

    printf STDERR "\nThis will be L2_rel_L1L2:\n";
    $derived_transform->{"L2_rel_L1L2"}=get_transform($data,$l2p,$forward);
    
    printf STDERR "\nThis will be L1L2_rel_L2:\n";
    $derived_transform->{"L1L2_rel_L2"}=get_transform($data,reverse_matched_pairs($l2p),$backward);
}


printf STDERR "\nThis will be L2_rel_L1L2:\n";
$derived_transform->{"L2_rel_L1L2"}=get_transform($data,$l2p,$forward);
# $derived_transform->{"L2_rel_L1L2"}=get_transform($data,$l2p,$backward);
# now try the same thing for L2 within the L1L2 assembly:
my @L2_in_L1L2_labels=@L2labels;
# transform labels in place
map {$_ =~ s/BS2535589/L1L2assy_BS2535589/} @L2_in_L1L2_labels;
# populate the coordinates of the new list according to the transform
populate_transformed_coords($data,[@L2_in_L1L2_labels],
			    $derived_transform->{"L2_rel_L1L2"},[@L2labels]);
# print out the rms 3-d variance between the transformed SMR group wrt the original SMR group.
printf STDERR "check: N = %d; y = %f\n",scalar(@L2_in_L1L2_labels),sqrt(constellation_var($data,mux_into_pairs([@L2_in_L1L2_labels],[@L2labels])));
# printf STDERR "%s\n",join("\n",@L2_in_L1L2_labels);

# how to keep track of the OCFs? as a separate set of "virtual" SMRs that form a triad? 4 would be required..
# origin, x-uvec, y-uvec, z-uvec = 4.
# triads are large here for the benefit of the fitting routine, which seems to run into a problem if triads are unity in dimension.
my $triad = {("origin" => {"X" =>    0, "Y" =>    0, "Z" =>    0},
	      "x_uvec" => {"X" => 1000, "Y" =>    0, "Z" =>    0},
	      "y_uvec" => {"X" =>    0, "Y" => 1000, "Z" =>    0},
	      "z_uvec" => {"X" =>    0, "Y" =>    0, "Z" => 1000})};
my $tran;

my @triad_labels=sort keys %{$triad};
# practice: output triads for each of the OCFs (L{1,2}S{1,2}):

# L1S1:
my @L1S1_OCF_triad_labels=@triad_labels;
map {$_ =~ s/^/L1S1_OCF_/} @L1S1_OCF_triad_labels;
printf STDERR "trans: %s\n",join(' ',@L1S1_OCF_triad_labels);
$tran={("dir" => 1,  "tz"  => 0)};
populate_transformed_coords($triad,[@L1S1_OCF_triad_labels],[$tran,$derived_transform->{"L1_rel_L1L2"}],[@triad_labels]);

# L1S2:
my @L1S2_OCF_triad_labels=@triad_labels;
map {$_ =~ s/^/L1S2_OCF_/} @L1S2_OCF_triad_labels;
printf STDERR "trans: %s\n",join(' ',@L1S2_OCF_triad_labels);
$tran={("dir" => 1,  "tz"  => +82.31)};
populate_transformed_coords($triad,[@L1S2_OCF_triad_labels],[$tran,$derived_transform->{"L1_rel_L1L2"}],[@triad_labels]);

# L2S1:
my @L2S1_OCF_triad_labels=@triad_labels;
map {$_ =~ s/^/L2S1_OCF_/} @L2S1_OCF_triad_labels;
printf STDERR "trans: %s\n",join(' ',@L2S1_OCF_triad_labels);
$tran={("dir" => 1,  "tz"  => -30.05)};
populate_transformed_coords($triad,[@L2S1_OCF_triad_labels],[$tran,$derived_transform->{"L2_rel_L1L2"}],[@triad_labels]);
# L2S2:
my @L2S2_OCF_triad_labels=@triad_labels;
map {$_ =~ s/^/L2S2_OCF_/} @L2S2_OCF_triad_labels;
printf STDERR "trans: %s\n",join(' ',@L2S2_OCF_triad_labels);
$tran={("dir" => 1,  "tz"  => 0)};
populate_transformed_coords($triad,[@L2S2_OCF_triad_labels],[$tran,$derived_transform->{"L2_rel_L1L2"}],[@triad_labels]);

# The coordinate system that the L{1,2}S{1,2} triads are reported relative to is the coordinate system of the TWE.
# so - in order to get the L1-L2 performance as described in the delivery package, the on-axis beam is given in this coordinate system.
# L1 system has r{x,y,z} = (-0.2, -0.2,+0.4  )microrad, t{x,y,z}=(0   ,-0.2 ,+0.2 )micron .. wrt OCF system
# L2 system has r{x,y,z} = (118.2,-3.8,+456.2)microrad, t{x,y,z}=(70.1,228.3,XX.XX)micron .. wrt OCF system

# OK!

# Move on to L3 measurements and figure where the OCFs will be when under operation.
# list of L3 SMRs relative to DATUM A is contained in: @L3labels_meas
# In 5351-110-int1-L-lens-assembly.pdf, page 1, note showing the apex measurements say:
# axial distance between L3S2 APEX and A:
# theoritical(sic): 2.975 mm
# after assembly: 3.027 mm
# under working operation: 3.27 mm *** (<- use this one)
# mechanical decenter of the L3 lens relative to the barrel as described in:
#    '5351-000-REP-141-01 - LSST L3 Test & measurement report.pdf'
# suggests that the L3 apex is translated in the XY plane and is located at
# 95 micron in X_meca
# 162 micron in Y_meca
# where there is a 3 degree clocking angle between the mechanical and CCS axes:
# \hat{X}_{meca} = cos(3deg)\hat{X}_{ccs} + sin(3deg)\hat{Y}_{ccs}
# \hat{y}_{meca} = cos(3deg)\hat{y}_{ccs} - sin(3deg)\hat{X}_{ccs}
# so expressed in CCS, the center of L3 is located at: 86.4e-3\hat{X}_{ccs} + 166.7e-3\hat{Y}_{ccs}.

my $L3_mech_decen={("dir" => -1,    "tx"  => 0.0864,    "ty"  => 0.1667)};

# L3S2:
my @L3S2_OCF_triad_labels=@triad_labels;
map {$_ =~ s/^/L3S2_OCF_/} @L3S2_OCF_triad_labels;

printf STDERR "trans: %s\n",join(' ',@L3S2_OCF_triad_labels);
$tran={("dir" => 1,  "tz"  => sum(-82.310,-412.642,-30.05,-346.237,0,-17.90,-54.1,-60))};
populate_transformed_coords($triad,[@L2S2_OCF_triad_labels],[$tran],[@triad_labels]);


# do this twice: once offset from L2S2 apex (assuming this follows the L2_CELL)
# and once offset from L1S1 apex..

printf STDERR "Two different representations for the L3 interface:\n\t'a' => transformed from L2S2;\n\t'b' => transformed from L1S1.\n";
my @L3IF_a_triad_labels=@triad_labels;
map {$_ =~ s/^/L3IF_a_/} @L3IF_a_triad_labels;
printf STDERR "trans: %s\n",join(' ',@L3IF_a_triad_labels);
$tran={("dir" => 1,  "tz"  => +3.27 + sum(-346.237,0,-17.90,-54.1,-60)-3.27)};

populate_transformed_coords($triad,[@L3IF_a_triad_labels],$tran,[@triad_labels]);
# L2S2_rel_L3S2 is the sum of thicknesses of layers 14..18 inclusive (e.g., config 3)
my $L2S2_rel_L3IF={%{$tran}};
printf STDERR "tz = %f\n",$L2S2_rel_L3IF->{"tz"};

$derived_transform->{"L3_a"}=get_transform($triad,mux_into_pairs([@L2S2_OCF_triad_labels],
								 [@L3IF_a_triad_labels]),$forward);
populate_transformed_coords($triad,[@L3IF_a_triad_labels],$derived_transform->{"L3_a"},[@triad_labels]);

my @L3IF_b_triad_labels=@triad_labels;
map {$_ =~ s/^/L3IF_b_/} @L3IF_b_triad_labels;
printf STDERR "trans: %s\n",join(' ',@L3IF_b_triad_labels);
# interpretation of tz: sum of +z sag of L3 apex wrt L3IF,
# (-z) thickness of L1, 
# (-z) separation between L1S2 & L2S1,
# (-z) thickness of L2,
# (-z) separation between L2S2 & FS1 (r band),
# (-z) thickness of F(r band),
# (-z) separation between FS2(r band) & L3S1,
# (-z) thickness of L3:

$tran={("dir" => 1,  "tz"  => +3.27 + sum(-82.310,-412.642,-30.05,-346.237,0,-17.90,-54.1,-60) -3.27)};

populate_transformed_coords($triad,[@L3IF_b_triad_labels],$tran,[@triad_labels]);
# L1S1_rel_L3S2 is the sum of thicknesses of layers 11..18 inclusive (e.g., config 3)
my $L1S1_rel_L3IF={%{$tran}};

printf STDERR "tz = %f\n",$L1S1_rel_L3IF->{"tz"};

$derived_transform->{"L3_b"}=get_transform($triad,mux_into_pairs([@L1S1_OCF_triad_labels],
								 [@L3IF_b_triad_labels]),$forward);
populate_transformed_coords($triad,[@L3IF_b_triad_labels],$derived_transform->{"L3_b"},[@triad_labels]);

printf STDERR "\nfinally the transformation between L3IF aligned with L2S2 vs. L3IF aligned with L1S1 is:\n";

get_transform($triad,mux_into_pairs([@L3IF_a_triad_labels],[@L3IF_b_triad_labels]),$forward);

# at this point, L3IF_a is with L3IF aligned to L2S2; L3IF_b is with L3IF aligned to L1S1.

# now try transforming the L3 SMRs (all measured WRT DATUM A [the L3IF]) into the coordinate system
# of the L1-L2 bundle:

my @L3labels_wrt_L1L2=@L3labels_meas;
map {$_ =~ s/^/L3_/;$_ =~ s/A3/L1L2assy/} @L3labels_wrt_L1L2;
# try this out: tip L3 (say, according to an intentional rotation, for a few angles.

my @L3_labels_wrt_L1L2_array=();

# this is where the L3 & Cryostat may be transformed (first rotated about the L3IF origin, then offset)
# and all SMR positions follow.
my $intentional_transform={};
$intentional_transform={("dir"=> 1,
			 "rx" => 0,		 "ry" => 0.1,		 "rz" => 0.0,
			 "tx" => 300.0,		 "ty" => 0.0,		 "tz" => 0.0)};
# null transform
$intentional_transform={("dir"=> 1,
			 "rx" => 0.0,		 "ry" => 0.0,		 "rz" => 0.0,
			 "tx" => 0.0,		 "ty" => 0.0,		 "tz" => 0.0)};

my @rot_L3_set=@L3labels_wrt_L1L2;
# map {$_ =~ s/$/_rot_$try_theta/} @rot_L3_set;


populate_transformed_coords($data,[@rot_L3_set],$intentional_transform,[@L3labels_meas]);
# and then place it into the L1L2 frame:
populate_transformed_coords($data,[@rot_L3_set],$derived_transform->{"L3_b"},[@rot_L3_set]);
push(@L3_labels_wrt_L1L2_array,@rot_L3_set);

# I thought I came across positions for the "visible" SMRs for L1, but maybe not. Need to look thru notes,
# can't find any code or data files that I thought I generated. for now we'll assume they don't exist, or we'll
# measure their positions at the same time as when alignment is performed between L1L2 & L3.

# what.qdp will contain the full set of SMR locations in the same coordinate system 
# as where the subset (L1L2) were measured at AOS (including those inaccessible) but all
# are now transformed (according the derived one using separate sets of measurements). So this set is
# the "best guess" of where the precisely measured SMRs now sit in the assembly.

save_measurements_file($data,[@L1_in_L1L2_labels,@L2_in_L1L2_labels,@L3_labels_wrt_L1L2_array],"what.qdp",0);

if ($output_pair_separations) {
    printf "\n";
    printf "2-SMR SEPARATION DISTANCES (%s):\n",join(',',"L2_in_L1L2_labels","L3_labels_wrt_L1L2_array");
    printf "\n";
    my $smr_groups=[
	[ map {(/Cell/)?():($_)} sort {(split("_",$a))[1] <=> (split("_",$b))[1]} @L2_in_L1L2_labels ],
	[ sort {(split("_",$a))[2] <=> (split("_",$b))[2]} @L3_labels_wrt_L1L2_array]
	];

    foreach my $L2_label (@{$smr_groups->[0]}) {
	foreach my $L3_label (@{$smr_groups->[1]}) {
	    printf "%s %s: %f\n",$L2_label,$L3_label,scalar_separation($data->{$L2_label},$data->{$L3_label});
	}
	printf "\n";
    }
}
# no more transformations needed, the L3 SMRs are already reported relative to DATUM A.

# where should DATUM A sit relative to L1S1 ??
# Well, it actually needs to sit somewhere relative to L2 because that has more power.. <<-- NOT TRUE!
# get the answer this way: set up a triad in the L3S2_OCF system that should coincide with L2S2_OCF.
# extract the transformation, then apply that transformation to the L3 SMRs to get a full representation of L3 in the L1L2 system.

# This derived transform (L3S2_rel_L2S2) expresses how the L3S2 apex would be transformed into
# the L1S1 coordinate system where it has been aligned precisely relative to L2S2:
# pars:
# (dir	=>	1)
# (rx	=>	0.000118237)
# (ry	=>	-3.81859e-06)
# (rz	=>	0.000456171)
# (tx	=>	0.068256)
# (ty	=>	0.171714)
# (tz	=>	1003.16)
# var(y) = 0.000000; rms deviation 0.000000
#
# BUT adding up the spacings in LSST_Ver_3.11_Baseline_Design
# the tz value (sum of thicknesses for layers 11..18) works out to 1003.239, which is 79 microns different.
# Need to check other versions of the design Zemax files.
# 1003.239 for baseline_design (79 microns different)
# 1003.307 for nonaxisymmetric_mod
# 1003.179 for nonaxisymmetrics_removed


# $data->{"jbldarf"}=$L2S2_rel_L3S2;

# print the whole label: not including the last argument causes a truncation of the labels delimited by "_".
save_measurements_file($triad,[map {($_ =~ /origin/)?($_):()} @L1S1_OCF_triad_labels,
			       map {($_ =~ /origin/)?($_):()} @L1S2_OCF_triad_labels,
			       map {($_ =~ /origin/)?($_):()} @L2S1_OCF_triad_labels,
			       map {($_ =~ /origin/)?($_):()} @L2S2_OCF_triad_labels,
			       map {($_ =~ /origin/)?($_):()} @L3IF_a_triad_labels,
			       map {($_ =~ /origin/)?($_):()} @L3IF_b_triad_labels],"whut.qdp",0);  

exit;
# now get the transformation needed to express the L1L2 composite in the (BALL) L2 frame

my $l1l2p = match_pairs([@L2labels],[@L1L2labels],["SMR\\d+","L\\d_Cell"]);

$derived_transform->{"L2subset_rel_L2Ball"}=get_transform($data,$l1l2p);


# make up a pairs list out of @L2labels and @L2_in_L1L2_labels to derive a transformation
# that will align the L1L2(AOS) set into the L2(Ball) L2S2 apex system:

my $l1l2p_l2s2 = mux_into_pairs([@L2labels],[@L2_in_L1L2_labels]);

# turns out the following transform is 100% identical to the one above, "L2subset_rel_L2Ball":
# but we'll use it anyway
$derived_transform->{"L2subset_rel_L2Ball2"}=get_transform($data,$l1l2p_l2s2);

# apply the populated data above (@L1_in_L1L2_labels & @L2_in_L1L2_labels) into this new frame

my @L1_in_L3if_labels=@L1_in_L1L2_labels;
map {$_ =~ s/L1L2assy/L3if/} @L1_in_L3if_labels;
my @L2_in_L3if_labels=@L2_in_L1L2_labels;
map {$_ =~ s/L1L2assy/L3if/} @L2_in_L3if_labels;

my @L1rot_in_L3if_labels=@L1_in_L1L2_labels;
map {$_ =~ s/L1L2assyrot/L3if/} @L1rot_in_L3if_labels;
my @L2rot_in_L3if_labels=@L2_in_L1L2_labels;
map {$_ =~ s/L1L2assyrot/L3if/} @L2rot_in_L3if_labels;

$derived_transform->{"L2subset_rel_L3if"}=$derived_transform->{"L2subset_rel_L2Ball"};
# this is where the additional offsets & rotations (should be included prior to computing SMR locations for L1 & L2 relative to L3 interface)
$derived_transform->{"L2subset_rel_L3if"}->{"tz"} += -478.580;


populate_transformed_coords($data,[@L1_in_L3if_labels,@L2_in_L3if_labels],
			    $derived_transform->{"L2subset_rel_L3if"},
			    [@L1_in_L1L2_labels,@L2_in_L1L2_labels]);

if (1) {
    $derived_transform->{"L2subset_rel_L3if"}->{"ry"} += 10*atan2(1,1)/45.0;

    populate_transformed_coords($data,[@L1rot_in_L3if_labels,@L2rot_in_L3if_labels],
				$derived_transform->{"L2subset_rel_L3if"},
				[@L1_in_L1L2_labels,@L2_in_L1L2_labels]);

    save_measurements_file($data,[@L1_in_L3if_labels,@L2_in_L3if_labels,
				  @L1rot_in_L3if_labels,@L2rot_in_L3if_labels,
				  @L3labels_meas],"L1L2L3SMRs_rel_L3if_datumA_rot.qdp");
    exit;
}

# and print out the transformed coords
save_measurements_file($data,[@L1_in_L3if_labels,@L2_in_L3if_labels,@L3labels_meas],"L1L2L3SMRs_rel_L3if_datumA.qdp");
# report transformations for each optical element according to this layout
# make up pairs for L1, L2..
my $l1_l3if = mux_into_pairs([@L1_in_L3if_labels],[@L1labels]);
my $l2_l3if = mux_into_pairs([@L2_in_L3if_labels],[@L2labels]);
my $dtf;
my $kys=["tx","ty","tz","rx","ry","rz"];

$dtf=get_transform($data,$l1_l3if);
printf "derived L1 tf wrt L3if: (%s)=(%s)\n",join(',',@{$kys}),join(',',@{$dtf}{@{$kys}});
$dtf=get_transform($data,$l2_l3if);
printf "derived L2 tf wrt L3if: (%s)=(%s)\n",join(',',@{$kys}),join(',',@{$dtf}{@{$kys}});

exit;
my @L1labels_nom = @L1labels;
my $got_transform;
map {$_ =~ s/BS2560710/nom_BS2560710/} @L1labels_nom;
populate_transformed_coords($data,[@L1labels_nom],$got_transform,[@L1labels]);


$got_transform=get_transform($data,$l2p);
my @L2labels_nom = @L2labels;
map {$_ =~ s/BS2535589/nom_BS2535589/} @L2labels_nom;
populate_transformed_coords($data,[@L2labels_nom],$got_transform,[@L2labels]);

# skip over the snippet below that just compares "as built"(meas) and "designed"(theo) pattern
# of SMRs on the L3 support ring, expressed in the interface coordinate system.

if (0) {
    printf STDERR "matches for L3 labels:\n";
    my $l3p = match_pairs([@L3labels],[@L3labels_meas],["\\d+_REF"]);

    printf "convar: %f\n",sqrt(constellation_var($data,$l3p,{"rx" => 0,"ry" => 0,"rz" => 0,"tx" => 0,"tz" => 0}));
    $got_transform=get_transform($data,$l3p);
}

# start out with some placeholder target values according to LSST optical design v3.3..
# L3 interface plane axial position wrt L3 apices: refer to 5351-110-int1-L-lens-assembly.pdf
# Figure of surfaces: 
# S2 radius is SR 13671.5 +/- 13mm. Expected under load: R 13360+/-3mm
# S1 radius at center is 3152.200 +/- 2mm, CA 722mm. Expected under load: R 3169+/-3mm
# S2 apex lies closer to FP than datum A: theoritical 2.975 after assembly 3.027 mm.
#    under working operation, 3.270 mm
# thickness of PEEK ring (part [4], "Axial Stop Ring") is a spacer between L3 & L3support
# and so is included in the 3.027 mm value of "S2 Apex to A" value. 
#    (PEEK shim theoretical thickness is 5.77mm)
# Data from optical design v.3.3 [check on these and get proper version number]
# FP dz   = -28.50mm
# L3 dz21 = -60mm
# L3 dz   = -53.300mm
# F/r dz21= -17.9mm
# F/r dz  = -349.58mm
# L2 dz21 = -30mm
# L2 dz   = -412.642mm
# L1 dz21 = -82.23mm
# L1 dz   = -3630.5mm # not relevant, this is distance from M3 apex.
#
# argument construction in lsst_mirrors.perl gives:
# asphere  -T -t 0 0 0 -- -R -19835 -k -1.215 -z 0 -ri 2558 -ro 4180 -A6 1.381e-09 
# | asphere  -T -t 0 0 -6156.2007 -- -R -6788 -k -0.222 -z 0 -ri 900 -ro 1710 -A6 -1.274e-05 -A8 -9.68e-07 
# | asphere  -T -t 0 0 233.7993 -- -R -8344.5 -k 0.155 -z 0 -ri 550 -ro 2508 -A6 -4.5e-07 -A8 -8.15e-09 
# | tran_ray -T -t 0 0 -11225 -- 
# | asphere -L  -T -t 0 0 7827.5267 -- -R1 -2824 -z1 0 -h1 20 -t1 20 -p1 1 -R2 -5021 -z2 -82.23 -ri 0 -ro 775 -p2 1 -h2 20 -t2 20 
# | asphere -L  -T -t 0 0 7332.6547 -- -R1 1e+34 -k1 0 -z1 0 -p1 1 -h1 20 -t1 20 -R2 -2529 -k2 -1.57 -z2 -30 -ri 0 -ro 551 -A6 0.001656 -p2 1 -h2 20 -t2 20 
# | asphere -L  -T -t 0 0 6953.0747 -- -R1 -5632 -z1 0 -p1 1 -h1 20 -t1 20 -R2 -5606 -z2 -17.9 -ri 0 -ro 375 -p2 1 -h2 20 -t2 20 -I 5629.995 7075.95 
# | asphere -L  -T -t 0 0 6884.0747 -- -R1 -3169 -k1 -0.962 -z1 0 -p1 1 -h1 20 -t1 20 -R2 13360 -z2 -60 -ri 0 -ro 361 -t2 20 -h2 20 -p2 1 
# | tran_ray -T -t 0 0 6795.574700 -- 
# | plane_collapse -c 0 0 1 0
#
# Separation between L2S2 and datum A (7332.6547 -30) - (6884.0747 -60) - 3.270mm =  475.310mm
# (so that under operation L2S2 to L3S2 is a distance 478.580mm apart)

# explore use of a hexapod adjustment algorithm. Specify a connection pattern and a standoff
# distance between two plates. Then derive the deltaL values needed to achieve a desired 
# inplane_rot-tip-tilt-decenterX-decenterY-piston adjustment *at a specific position*. 
# in the lab, two iterations will probably be necessary especially since the junction 
# points of the hexapod in 3space aren't truly known, they're estimated from this ideal model.

my $fixed_hp_radius=500.0;
my $unfixed_hp_radius=700.0;
my $fixed_hp_zcoord=-200.0;   # absolute coord
my $unfixed_hp_zcoord=+100.0; # absolute coord
my $hp_adj_pt={"X"=>0,
	       "Y"=>0,
	       "Z"=>+500.0};

my $fixed_hp_azlist=[(0-15),(0+15),
		     (120-15),(120+15),
		     (240-15),(240+15)];
my $unfixed_hp_azlist=[ map {60+$_} @{$fixed_hp_azlist} ];

printf STDERR "pre-roll fixed  : %s\n",join(',',@{$fixed_hp_azlist});
printf STDERR "pre-roll unfixed: %s\n",join(',',@{$unfixed_hp_azlist});
# roll the first element in $unfixed_hp_azlist
{
    my $elem=splice(@{$unfixed_hp_azlist},$#{$unfixed_hp_azlist},1); # pop
    splice(@{$unfixed_hp_azlist},0,0,$elem); # unshift
}
printf STDERR "post-roll fixed  : %s\n",join(',',@{$fixed_hp_azlist});
printf STDERR "post-roll unfixed: %s\n",join(',',@{$unfixed_hp_azlist});

my $struts=[];

foreach my $ix (0..$#{$fixed_hp_azlist}) {
    # populate "fixed" and "unfixed_nom" fields of each element of $struts
    $struts->[$ix]={"fixed" => {"X" => $fixed_hp_radius*cos($fixed_hp_azlist->[$ix]*$deg),
				"Y" => $fixed_hp_radius*sin($fixed_hp_azlist->[$ix]*$deg),
				"Z" => $fixed_hp_zcoord},
		    "unfixed_nom" => {"X" => $unfixed_hp_radius*cos($unfixed_hp_azlist->[$ix]*$deg),
				      "Y" => $unfixed_hp_radius*sin($unfixed_hp_azlist->[$ix]*$deg),
				      "Z" => $unfixed_hp_zcoord}};
#    printf STDERR "fixed  : %s\n",show_vec($struts->[$ix]->{"fixed"});
#    printf STDERR "unfixed: %s\n",show_vec($struts->[$ix]->{"unfixed_nom"});
    $struts->[$ix]->{"length_nom"}=scalar_separation($struts->[$ix]->{"fixed"},$struts->[$ix]->{"unfixed_nom"});
#    printf STDERR "strut length (%d) = %f\n",$ix,$struts->[$ix]->{"length_nom"};
}

# now compute each "unfixed_adjust" point based on transforming the points evaluated relative to the $hp_adj_pt
foreach my $ix (0..$#{$unfixed_hp_azlist}) {
    $struts->[$ix]->{"unfixed_adjust"}=vec_separation($struts->[$ix]->{"unfixed_nom"},$hp_adj_pt);
    my $tf={"rx" => 0,"ry" => 0,"rz" => 0*$deg,"tx" => 0,"ty" => 0,"tz" => 1};
    printf STDERR ("unfixed_nom: %s (dl=%f)\n",
		   show_vec($struts->[$ix]->{"unfixed_adjust"}),
		   scalar_separation(vec_separation($struts->[$ix]->{"unfixed_adjust"},flip($hp_adj_pt)),$struts->[$ix]->{"fixed"})-
		   $struts->[$ix]->{"length_nom"});

    $struts->[$ix]->{"unfixed_adjust"}=transform($struts->[$ix]->{"unfixed_adjust"},$tf);

    $struts->[$ix]->{"strut_adj"}=
	scalar_separation(vec_separation($struts->[$ix]->{"unfixed_adjust"},
					 flip($hp_adj_pt)),
			  $struts->[$ix]->{"fixed"})-
			  $struts->[$ix]->{"length_nom"};
    
    printf STDERR ("unfixed_adj: %s (dl=%f)\n",show_vec($struts->[$ix]->{"unfixed_adjust"}),$struts->[$ix]->{"strut_adj"});
#    printf STDERR "unfixed_adj: %s (dl=%f)\n",show_vec($struts->[$ix]->{"unfixed_adjust"});
}

# now move on to place the L1/L2 assembly aligned with L3/cryostat interface plane. (tip L3/Cryostat interface and figure out how to follow L3 with L1/L2
# use only the THALES report 5351-IR 0005-0011012639 to define the L3/cryostat interface plane
# The points list is represented as @{$data}{@L3labels_meas} with (X,Y,Z) aligned with CCS system.
# The required shift to the @{$data}{@L1L2labels} needed to establish design spacing 
# between loaded L3S2 apex and L2S2 apex is: delta Z = -478.580 mm (see above calculation)
# make a new slice kept in $data:
if (0) {
    foreach my $entry ( @L1L2labels ) {
	printf STDERR "entry $entry\n";
    }
}


exit;

# get regression for each set of SMRs for the nearest planar fit
my $face_smrs=Statistics::Regression->new("face SMR model",["z0","dz/dx","dz/dy"]);
my $radial_smrs=Statistics::Regression->new("radial SMR model",["z0","dz/dx","dz/dy"]);
my $smrvec=[$face_smrs,$radial_smrs];
my $smrgrp=[];


foreach my $smr_entry ( @labels ) {
    next if ($smr_entry !~ /REF/);
    next if ($smr_entry =~ /THEO/);
    my $ix=(split("_",$smr_entry))[1];
    my $smix=int(($ix-1)/12);
    $smrgrp->[$smix]=[] if (!defined($smrgrp->[$smix]));
    push(@{$smrgrp->[$smix]},$smr_entry);
    $smrvec->[$smix]->include($data->{$smr_entry}->{"Z"},[1.0,@{$data->{$smr_entry}}{"X","Y"}]);
}

# get back the regressions
foreach my $smix (0,1) {
    $smrvec->[$smix]->print();
    my $theta=$smrvec->[$smix]->theta();
    printf STDERR "model: %s\n",join(' ',@{$theta});
    printf STDERR "avg. deviation: %g\n",sqrt($smrvec->[$smix]->sigmasq()/$smrvec->[$smix]->n());
}

open(E,">","L3SMR_measurement_analysis.qdp") || die;
printf E ("lag g1 X[mm]\n");
printf E ("lag g2 Y[mm]\n");
printf E ("lag g3 Z[mm]\n");
printf E ("lag g4 X\\dmeas\\u-X\\dtheo\\u[mm]\n");
printf E ("lag g5 Y\\dmeas\\u-Y\\dtheo\\u[mm]\n");
printf E ("lag g6 Z\\dmeas\\u-Z\\dtheo\\u[mm]\n");
printf E ("lag g7 Z\\dmeas\\u-Z\\dmodel\\u[mm]\n");

foreach my $smix (0,1) {
    foreach my $smr_entry ( @{$smrgrp->[$smix]} ) {
	next if ($smr_entry =~ /THEO/);

	my $theta=$smrvec->[$smix]->theta();
	my $zmodel = $theta->[0];
	$zmodel   += $theta->[1]*$data->{$smr_entry}->{"X"};
	$zmodel   += $theta->[2]*$data->{$smr_entry}->{"Y"};

	my $smr_theo=$smr_entry;
	$smr_theo =~ s/MEAS/THEO/;
	printf E ("%f %f %f %f %f %f %f\n",
		  @{$data->{$smr_entry}}{"X","Y","Z"},
		  $data->{$smr_entry}->{"X"}-$data->{$smr_theo}->{"X"},
		  $data->{$smr_entry}->{"Y"}-$data->{$smr_theo}->{"Y"},
		  $data->{$smr_entry}->{"Z"}-$data->{$smr_theo}->{"Z"},
		  $data->{$smr_entry}->{"Z"}-$zmodel,
#		  $data->{$smr_entry}->{"Z"},
#		  $zmodel
	    );
    }
    printf E "no no no no no no no\n";
}
close(E);

exit;

sub save_measurements_file {
    my ($data,$entrylist,$outputfile,$startix)=@_;
    my $labix=10;
    printf STDERR " .. %s .. ",$outputfile;
    open(Q,">",$outputfile) || die;
    foreach my $smr_entry ( @{$entrylist} ) {
	my @lab=split("_",$smr_entry);
	my $lab;
	$startix=(defined($startix))?$startix:2;
	$lab=join("_",@lab[$startix..$#lab]);
	printf Q "%f %f %f\n",@{$data->{$smr_entry}}{"X","Y","Z"};
	printf Q "lab %d csi 0.7 pos %f %f to %f %f jus lef \"%s\"\n",$labix,@{$data->{$smr_entry}}{"X","Y"},$data->{$smr_entry}->{"X"}+20,$data->{$smr_entry}->{"Y"}+20,$lab;
	$labix++;
    }
    close(Q);
}

sub transform {
    # 
    # this routine has been modified to accept either a single transformation directive hash or an array of them.
    # transformations are performed in the order they are received: $tf_list->[0], $tf_list->[1] .. etc.
    # 
    my ($p,$tf_list)=@_;
    my $tfp={%{$p}};
    my ($cs,$sn);
    # do a 3-2-1 transformation followed by a straight translation

    if (ref($tf_list) eq "HASH") {
	# wrap the hash into a list and use below as a list would be.
	$tf_list = [$tf_list];
    }

    foreach my $tf (@{$tf_list}) {
	
	die "need to specify transformation direction!!\n" if (!defined($tf->{"dir"}));
	if ($tf->{"dir"}==1) { # forward transform, as if taking the template, rotating and then translating it to match the data
	    if (defined($tf->{"rz"})) {
		($cs,$sn)=(cos($tf->{"rz"}),sin($tf->{"rz"})); 
		$tfp={"X" => $cs*$tfp->{"X"}-$sn*$tfp->{"Y"},
			  "Y" => $cs*$tfp->{"Y"}+$sn*$tfp->{"X"}, 
			  "Z" => $tfp->{"Z"}};
	    }
	    
	    if (defined($tf->{"ry"})) {
		($cs,$sn)=(cos($tf->{"ry"}),sin($tf->{"ry"})); 
		$tfp={"Z" => $cs*$tfp->{"Z"}-$sn*$tfp->{"X"},
			  "X" => $cs*$tfp->{"X"}+$sn*$tfp->{"Z"}, 
			  "Y" => $tfp->{"Y"}};
	    }
	    
	    if (defined($tf->{"rx"})) {
		($cs,$sn)=(cos($tf->{"rx"}),sin($tf->{"rx"}));
		$tfp={"Y" => $cs*$tfp->{"Y"}-$sn*$tfp->{"Z"},
			  "Z" => $cs*$tfp->{"Z"}+$sn*$tfp->{"Y"}, 
			  "X" => $tfp->{"X"}};
	    }
	    
	    $tfp->{"X"} += $tf->{"tx"} if (defined($tf->{"tx"}));
	    $tfp->{"Y"} += $tf->{"ty"} if (defined($tf->{"ty"}));
	    $tfp->{"Z"} += $tf->{"tz"} if (defined($tf->{"tz"}));
	} elsif ($tf->{"dir"}==-1) { # reverse transform, as if taking the data, reversing the translation, reversing the rotation to compare to the template.
	    # reverse ordering and flip..
	    $tfp->{"Z"} -= $tf->{"tz"} if (defined($tf->{"tz"}));
	    $tfp->{"Y"} -= $tf->{"ty"} if (defined($tf->{"ty"}));
	    $tfp->{"X"} -= $tf->{"tx"} if (defined($tf->{"tx"}));
	    
	    if (defined($tf->{"rx"})) {
		($cs,$sn)=(cos(-1*$tf->{"rx"}),sin(-1*$tf->{"rx"}));
		$tfp={"Y" => $cs*$tfp->{"Y"}-$sn*$tfp->{"Z"},
			  "Z" => $cs*$tfp->{"Z"}+$sn*$tfp->{"Y"}, 
			  "X" => $tfp->{"X"}};
	    }
	    if (defined($tf->{"ry"})) {
		($cs,$sn)=(cos(-1*$tf->{"ry"}),sin(-1*$tf->{"ry"})); 
		$tfp={"Z" => $cs*$tfp->{"Z"}-$sn*$tfp->{"X"},
			  "X" => $cs*$tfp->{"X"}+$sn*$tfp->{"Z"}, 
			  "Y" => $tfp->{"Y"}};
	    }
	    if (defined($tf->{"rz"})) {
		($cs,$sn)=(cos(-1*$tf->{"rz"}),sin(-1*$tf->{"rz"})); 
		$tfp={"X" => $cs*$tfp->{"X"}-$sn*$tfp->{"Y"},
			  "Y" => $cs*$tfp->{"Y"}+$sn*$tfp->{"X"}, 
			  "Z" => $tfp->{"Z"}};
	    }
	} else {
	    die "transformation direction missing.\n";
	}
    }
    return($tfp);
}

sub constellation_var { # computes the variance between two sets of coordinates when a transformation is applied to the second of each.
    my ($data,$pairs,$tf,$outfile)=@_;
    my $chisq=0;
    my $n_pair=0;
    my $var_by_point={};

    open(F,">",$outfile) if (defined($outfile));
    foreach my $pair (@{$pairs}) {
	my ($d1,$d2)=@{$data}{@{$pair}};
	# if $tf is defined, make a copy of $d1 and transform in place.
	$d2=transform($d2,$tf) if (defined($tf));
	my $d={};

	my $pair_str=join(',',@{$pair});
	$var_by_point->{$pair_str}=0;

	foreach my $ax ("X","Y","Z") {
	    $d->{$ax}=$d1->{$ax}-$d2->{$ax};
	    my $thisdev=pow($d->{$ax},2.0);
	    $var_by_point->{$pair_str} += $thisdev;
	    $chisq += $thisdev
	#	    printf STDERR "pair diff $pair->[0] ($ax): %f\n",$d1->{$ax}-$d2->{$ax};
	}
	if (defined($outfile)) {
	    printf F "win 2\nlab %d pos %f %f ma 13 msiz 2 \" \"\n",10+2*$n_pair,@{$d1}{"X","Y"};
	    printf F "win 3\nlab %d pos %f %f ma 13 msiz 2 \" \"\n",10+2*$n_pair+1,@{$d1}{"X","Z"};
	    printf F "%f %f %f\n",@{$d1}{"X","Y","Z"};
	    printf F ("%f %f %f\n",$d1->{"X"}-1000*$d->{"X"},$d1->{"Y"}-1000*$d->{"Y"},$d1->{"Z"}-1000*$d->{"Z"});
	    printf F "no no no\n";
	}
	$n_pair++;
    }
    # output the maximum deviant if $outfile is defined. otherwise do nothing because same routine is used in iterative fitting
    if (defined($outfile)) {
	# find the largest offender
	my $worst_offenders = [ reverse sort {$var_by_point->{$a} <=> $var_by_point->{$b}} keys %{$var_by_point} ];
	foreach my $ix (0..$#{$worst_offenders}) {
	    printf F "! worst offender %d (%s) = %f\n",$ix+1,$worst_offenders->[$ix],sqrt($var_by_point->{$worst_offenders->[$ix]});
	}
    }
    close(F);
    return($chisq/$n_pair);
}

sub reverse_matched_pairs {
    my ($list)=@_;
    # return an array containing 2-arrays with orders flipped
    return( [ map { [reverse @{$_}] } @{$list}] );
}

sub match_pairs {
    my ($entry1,$entry2,$mtch)=@_;
    my $matches=[];

    foreach my $ent1 (@{$entry1}) {

	my $nomatch=0;
	my $match={};
	foreach my $mtc (@{$mtch}) {
	    ($match->{$mtc})=($ent1 =~ /_($mtc)_/);
	    $nomatch=1 if (!defined($match->{$mtc}));
	}
	next if ($nomatch);

	foreach my $ent2 (@{$entry2}) {
	    $nomatch=0;
	    foreach my $mtc (@{$mtch}) {
		$nomatch=1 if ($ent2 !~ /_$match->{$mtc}_/);
	    }
	    next if ($nomatch);
	    push(@{$matches},[$ent1,$ent2]);
	}
    }
    foreach my $ix (0..$#{$matches}) {
	printf STDERR "matching (e.g.) %d: %s\n",$ix,join(',',@{$matches->[$ix]})if ($ix == 0);
    }
    return($matches); # array containing pairs of entries (as 2-arrays)
}

my $globals={};

sub get_transform {
    # this routine returns a transform hash (not an array of them) based on (genetically)
    # modifying the transform that's passed to it.
    my ($data,$pairs,$tf)=@_;

    my @tf_fields=("tx","ty","tz","rx","ry","rz","dir");
    my %tf_field_errs=("tx"=>1,	       "ty"=>1,	       "tz"=>1,
		       "rx"=>0.01,     "ry"=>0.01,     "rz"=>0.01,
		       "dir"=>0);
    
    my $guess=[(0) x (scalar(@tf_fields)-1), $tf->{"dir"} ];
    my $scale=[@tf_field_errs{@tf_fields}];
    $scale=[@{$tf->{"scale"}}{@tf_fields}] if (defined($tf->{"scale"}));
    $globals->{"data"}=$data;
    $globals->{"tf_fields"}=[@tf_fields];
    $globals->{"pairs"}=$pairs;
#    printf STDERR "number of pairs: %d\n",scalar(@{$pairs});
#    map {printf STDERR "\t%s -> %s\n",@{$_}[0,1]} @{$pairs};
    my ($p,$y)=MinimiseND($guess,$scale,\&match_points,1e-8,1e6);
    my $final={};
    @{$final}{@tf_fields}=@{$p};
    printf STDERR "pars:\n%s\n",join("\n",map {sprintf("(%s\t=>\t%g)",$_,$final->{$_})} (sort @tf_fields));
#    join(',',sort @tf_fields),join(',',@{$final}{sort @tf_fields});
    printf STDERR "N = %d; var(y) = %f; rms deviation %f\n",scalar(@{$pairs}),$y,sqrt($y);
    return($final);
}

sub match_points {
    my $tf={};
    @{$tf}{@{$globals->{"tf_fields"}}}=@_;
    return(constellation_var($globals->{"data"},$globals->{"pairs"},$tf));
}

sub show_vec {
    my ($v)=@_;
    return(sprintf("(X,Y,Z)=(%f,%f,%f)",@{$v}{"X","Y","Z"}));
}

sub vec_modulus {
    my ($v)=@_;
    my $mod=0;
    foreach my $ax ("X","Y","Z") {
	$mod+=pow($v->{$ax},2);
    }
    return($mod);
}

sub invert_transform {
    my ($tf_list)=@_;
    my $mod_tf_list;
    if (ref($tf_list) eq "HASH") {
	$mod_tf_list = {%{$tf_list}};
	$mod_tf_list->{"dir"} *= -1;
    } elsif (ref($tf_list) eq "ARRAY") {
	$mod_tf_list=[];
	while (my $tf=pop(@{$tf_list})) {
	    my $mod_tf={%{$tf}};
	    $mod_tf->{"dir"} *= -1;
	    push(@{$mod_tf_list},$mod_tf);
	}
    } else {
	printf STDERR "tf_list: $tf_list\n";
	die "what kind of reference is this??\n";
    }

    $mod_tf_list;
}

sub flip {
    my ($v)=@_;
    $v = {%{$v}}; 
    foreach my $ax ("X","Y","Z") {
	$v->{$ax} *= -1;
    }
    return($v);
}

sub vec_separation {
    my ($p1,$p2)=@_;
    my $d={};
    foreach my $ax ("X","Y","Z") {
	$d->{$ax}=($p1->{$ax}-$p2->{$ax});
    }    
    return($d);
}

sub scalar_separation {
    my ($p1,$p2)=@_;
    my $d=vec_separation($p1,$p2);
    return(sqrt(vec_modulus($d)));
}

sub populate_transformed_coords {
    my ($data,$tf_entries,$tf,$orig_entries)=@_;

    if (ref($tf) eq "HASH") {
	printf STDERR ("will transform coordinates using:\n%s\n",
		       join("\n",map {sprintf("\t(%s\t=>\t%g)",$_,$tf->{$_})} (sort keys %{$tf})));
    } else {
	my @out=();
	foreach my $tf0 (@{$tf}) {
	    push(@out,join("\n",map {sprintf("\t(%s\t=>\t%g)",$_,$tf0->{$_})} (sort keys %{$tf0})));
	}
	printf STDERR (join("\n----------------------------------------------\n","will transform coordinates using:",@out,"\n"));
    }
    
    my $null_tran={"rx" => 0,"ry" => 0,"rz" => 0,"tx" => 0,"ty" => 0,"tz" => 0,"dir" => 1};
    my $sumvar=0;
    foreach my $ix (0..$#{$orig_entries}) {
	$data->{$tf_entries->[$ix]} = transform($data->{$orig_entries->[$ix]},$tf);
	$sumvar+=sum(map {$_*$_} (map {$data->{$tf_entries->[$ix]}->{$_}-$data->{$orig_entries->[$ix]}->{$_}} ("X","Y","Z")));
    }
    $sumvar /= scalar(@{$orig_entries});
    printf STDERR ("rms d(x,y,z): %g\n",sqrt($sumvar));
    return;
}

sub mux_into_pairs {
    my ($set1,$set2)=@_;
    die "sets have different lengths!" if ($#{$set1} != $#{$set2});
    my $pairs=[];
    foreach my $ix (0..$#{$set1}) {
	push(@{$pairs},[$set1->[$ix],$set2->[$ix]]);
    }
    return($pairs);
}

sub sum {
    my @list=@_;
    my $sum=0;
    map {$sum += $_} @list;
    return($sum);
}
