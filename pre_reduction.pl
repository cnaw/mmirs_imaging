#!/usr/bin/perl -w
#
# Code to create a Task Control file for MMIRS reduction
# spectroscopy section very incomplete
# 
# Christopher N. A. Willmer 
# cnaw@as.arizona.edu
# version 0.0  2017-05-18
#
use Astro::FITS::CFITSIO qw( :longnames :shortnames :constants PerlyUnpacking );
use Carp ;
#
if($#ARGV <= 2) {
    print "pre_reduction.pl date field filter crosstalk two-point:\n\n";
    print "pre_reduction.pl 2018.0829 sw Y 37.5 two_point\nor\n";
    print "pre_reduction.pl 2018.0829 sw Y 37.5\nor\n";
    print "pre_reduction.pl 2018.0829 sw Y\n";
    exit(0);
} else {
    $date   = $ARGV[0];
    $field  = $ARGV[1];
    $filter = $ARGV[2];
    if($#ARGV == 3) {
	$crosstalk = $ARGV[3];
    } else {
	$crosstalk = 0;
    }
    if($#ARGV == 4) {
	$two_point = 1;
    } else {
	$two_point = 0;
    }
}
#$root_dir  = '/data/Local/cnaw/nep/';
#$programme = '2017b-UAO-S117';
# GOGREEN
#$root_dir = '/data1/Local/cnaw/uao_s113/' ;
# Y, H
$root_dir = '/data1/Local/cnaw/nep/' ;
$programme = '/' ;
$verbose = 0;

$field     = lc($field);

# dont forget to edit 
# /home/cnaw/idl/mmirs/pipeline/crosstalk_correction_32amp.pro
# line 17
# values can be 00 -25, -37.5, -50
# 
#$xtalk_dir = $field.'_y_37.5/';
if($crosstalk == 0) {
    $xtalk_dir = join('_',$field, lc($filter),'00/');
}
if($crosstalk == 25) {
    $xtalk_dir = join('_',$field, lc($filter),'25/');
}
if($crosstalk == 37.5) {
    $xtalk_dir = join('_',$field, lc($filter),'37.5/');
}
if($crosstalk == 50) {
    $xtalk_dir = join('_',$field, lc($filter),'50/');
}

#
# from here things should be fixed
#
#$template = $programme.'/'.$date.'/*.fits';
#$template  = $root_dir.$programme.$date.'/*$field*.fits';
$science   = $root_dir.$programme.$date.'/*'.uc($field).'*.fits';
$darks     = $root_dir.$programme.$date.'/dark*.fits';
$idl_batch = $root_dir.$date.'_reduce.idl';
$date_dir  = $date.'/';
#$year = substr($date,0,4);
@files   = `ls $science $darks | grep -v skycam`;
@science_type = ();
@flats   = ();
@arcs    = ();
print "@files\n";
print "$science\n$darks\npause";
<STDIN>;
#
#========================================================================
#
# Step 1.
# Read file headers and group according to image type 
# (object, arc, flat, dark) and exposure time (darktime)
#
for($l = 0 ; $l <= $#files ; $l++) {
    $files[$l] =~ s/\n//g;
    print "$files[$l]\n";
#
    $status    = 0;
    my($fptr) = Astro::FITS::CFITSIO::open_file($files[$l],Astro::FITS::CFITSIO::READONLY(),$status);
    check_status($status) or die;
#
# skip files without extensions (for imaging at least)
#
    ffgkey($fptr,'EXTEND',$extend,$comment,$status);
    if($status != 0) { next;}
#
# Move to first extension
# 
    ffmrhd($fptr,1,$hdutype,$status)    ;
    $status = 0;
    ffgkey($fptr,'EXPTIME',$exptime,$comment,$status);
    check_status($status) or die;

#    $status = 0;
#    ffgkey($fptr,'EXPTABLE',$exptable,$comment,$status);
#    check_status($status) or die;

    $status = 0;
    ffgkey($fptr,'DARKTIME',$darktime,$comment,$status);
    check_status($status) or die;
#
# recover filter value
#    
    $status = 0;
    ffgkey($fptr,'FILTER',$filter,$comment,$status);
    check_status($status) or die;
#
# find slit mask name or open in the case of imaging
#
    $status = 0;
    ffgkey($fptr,'APERTURE',$aperture,$comment,$status);
    check_status($status) or die;
#
# Grism or open
#
    $status = 0;
    ffgkey($fptr,'DISPERSE',$disperse,$comment,$status);
    check_status($status) or die;
#
    $fptr->close_file($status);
    if($verbose > 0) {print "$files[$l] $darktime $filter\n";}
#
# group files according to dark exposure time, storing file
# name in a hash array, using the dark time as key
# remove directory path from image name.
# 
    @junk = split('\/',$files[$l]);
    $files[$l] = $junk[$#junk];

    $darktime = sprintf("%8.3f",$darktime);
    if(exists($darks{$darktime})) {
	$darks{$darktime} = join(' ',$darks{$darktime}, $files[$l]);
    } else {
	$darks{$darktime} = $files[$l];
    }
    $info{$files[$l]} = join(' ',$exptime, $darktime, $filter, $aperture, $disperse);
#
# group files according to filter, which is the hash array key
#
    if($filter !~ m/dark/) {
	if(exists($filters{$filter})) {
	    $filters{$filter} = join(' ',$filters{$filter}, $files[$l]);
	} else {
	    $filters{$filter} = $files[$l];
	}
    }
# skip darks
    if($files[$l] =~ m/dark/) { next;}
# store arcs
    if($files[$l] =~ m/comp/) {
	push(@arcs, $files[$l]);
	next;
    }
# store flats
    if($files[$l] =~ m/flat/) {
	push(@flats, $files[$l]);
	next;
    }
#
# store science images, taking into account 
# if these are images or spectra
#
    push(@science, $files[$l]);
    push(@science_type, $aperture);
}
#
# print files sorted according to darktime
#
if($verbose > 0 ) {
    foreach $darktime (sort(keys(%darks))){
	@files = split(' ',$darks{$darktime});
	for($j = 0 ; $j <= $#files ; $j++){
	    print "$darktime $files[$j]\n";
	}
    }
    print "pause";
    <STDIN>;
}
#
#===========================================================
#
# Step 2.
# generate the task control files for science images. Each
# exposure should have its own file in the case of imaging
# and each pair of nods in the case of spectroscopy.
#
open(IDL, ">$idl_batch") || die "cannot open $idl_batch";
print "\nScience images\n\n";
for($file_in_list = 0 ; $file_in_list <= $#science ; $file_in_list++) {
    $file = $science[$file_in_list];
    if($science_type[$file_in_list] !~ m/open/) { next;}
    print "$file $science_type[$file_in_list]\n";
    $logfile = $file;
    $logfile =~ s/.fits/.txt/g;
    open(LOG,">$logfile") || die "cannot open $logfile";
    ($exptime, $darktime, $filter, $aperture, $disperse) = split(' ',$info{$file});

    $file =~ s/.fits//g;
#
    $keyword = 'SIMPLE';
    $value   = 'T';
    $comment = 'FITS-like header';
    $line = sprintf("%-8s= %-20s /%-20s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
    $keyword = 'LONGSTR';
    $value   = 'OGIP 1.0';
    $comment = 'The OGIP long string convention may be used';
    $line = sprintf("%-8s= %-20s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# raw data directory
#
    $keyword = 'RAW_DIR';
#    $value   = $root_dir.$programme.'/'.$date_dir;
    $value   = $root_dir.$date_dir;
    $comment = ' ';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# pre-processed data go here
#
    $keyword = 'R_DIR';
    $value   = $root_dir.'preproc/'.$date_dir.$xtalk_dir;
    $comment = ' ';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# reduced data go here
#
    $keyword = 'W_DIR';
    $value   = $root_dir.'reduced/'.$date_dir.$xtalk_dir;
    $comment = ' ';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# extension of raw data
#
    $keyword = 'RAWEXT';
    $value   = '.gz';
    $comment = ' ';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# instrument
#
    $keyword = 'INSTRUME';
    $value   = 'MMIRS';
    $comment = 'spectrograph name';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# grism
#
    $keyword = 'GRISM';
    $value   = $disperse;
    $comment = 'grism';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# filter
#
    $keyword = 'FILTER';
    $value   = $filter;
    $comment = 'filter';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# bright source
#
    $keyword = 'BRIGHT';
    $value   = 0 ;
    $comment = 'bright source 1=Yes 0=No';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#-----------------------------------------------------------------------
#
# Science image and associated darks
#
    $keyword = 'SCI';
    $value   = "'".$file."'" ;
    $comment = 'science frame (only one)';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# if this is a spectroscopy file, and the second nod and
# increase counter
    if($science_type[$file_in_list] !~ m/open/) {
	$file_in_list++;
	$file = $science[$file_in_list];
	$keyword = 'SCI2';
	$value   = "'".$file."'" ;
	$comment = 'science frame (only one)';
	$line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
	print LOG $line,"\n";
    }
#
# darks for science frame
#
    $darktime = sprintf("%8.3f",$darktime);
    @files = split(' ',$darks{$darktime});
    @darks = ();
    for($j = 0 ; $j <= $#files ; $j++) {
	if($files[$j] =~ m/dark/) {
	    push(@darks,$files[$j]);
	}
    }
    $list = "'".$darks[0];
    for($j = 1; $j <= $#darks; $j++) {
	$list = join(',',$list,$darks[$j]);
    }
    $list = $list."'";
    $list =~ s/.fits//g;
#
# need to check the total string length, and if longer than 68 characters
# create continuation line
#
    $keyword = 'DARKSCI';
    
    $len  = length($list);
    print "$keyword length is $len\n";
    if($len <= 69) {
	$value = $list ;
	$line = sprintf("%-8s= %-70s",$keyword, $value);
	print LOG $line,"\n";
    } else {
# break into two lines
	$value   =  substr($list,0,68)."&'";
	$line2   =  "'".substr($list,68);
	$comment = 'dark frames for science';
	$line = sprintf("%-8s= %-70s",$keyword, $value);
	print LOG $line,"\n";
#
	$keyword = 'CONTINUE';
	$line = sprintf("%-8s %-70s",$keyword, $line2);
	print LOG $line,"\n";
	
	$keyword = 'LONGSTRN';
	$value   = '';
	$line = sprintf("%-9s",$keyword);
	print LOG $line,"\n";
    }
#
#-----------------------------------------------------------------------
# flats
#
    $nflats = 0;
    $list = "'";
    for ($j = 0; $j <= $#flats ; $j++) {
	$flat_file = $flats[$j];
	@junk = split(' ',$info{$flat_file});
	if($junk[2] eq $filter) {
	    if($list eq "'") {
		$darktime = sprintf("%8.3f",$junk[1]);
		$list = $flat_file;
		$nflats++;
	    } else {
		$list = join(',',$list,$flat_file);
		$nflats++;
	    }
	}
    }
    if($nflats > 0) {
	$list = $list."'";
	$list =~ s/.fits//g;
	$keyword = 'FLAT';
	$value   =  $list;
	$comment = 'flat images';
	$line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
	print LOG $line,"\n";
#
# now the corresponding darks
#
	@darks = ();
	for($j = 0 ; $j <= $#files ; $j++) {
	    if($files[$j] =~ m/dark/) {
		push(@darks,$files[$j]);
	    }
	}
	$list = "'".$darks[0];
	for($j = 1; $j <= $#darks; $j++) {
	    $list = join(',',$list,$darks[$j]);
	}
	$list = $list."'";
	$list =~ s/.fits//g;
	$keyword = 'DARKFLAT';
	$len  = length($list);
	print "$keyword length is $len\n";
	if($len <= 69) {
	    $value = $list;
	    $line = sprintf("%-8s= %-70s",$keyword, $value);
	    print LOG $line,"\n";
	} else {
# break into two lines
	    $value   =  substr($list,0,68)."&'";
	    $line2   =  "'".substr($list,68);
	    $comment = 'dark frames for flats';
	    $line = sprintf("%-8s= %-70s",$keyword, $value);
	    print LOG $line,"\n";
#
	    $keyword = 'CONTINUE';
	    $line = sprintf("%-8s %-70s",$keyword, $line2);
	    print LOG $line,"\n";
	    
	    $keyword = 'LONGSTRN';
	    $value   = '';
	    $line = sprintf("%-9s",$keyword);
	    print LOG $line,"\n";
	}
    }
#-----------------------------------------------------------------------
# arcs
#
    $narcs = 0;
    $list = "'";
    for ($j = 0; $j <= $#arcs ; $j++) {
	$arc_file = $arcs[$j];
	@junk = split(' ',$info{$arc_file});
	if($junk[2] eq $filter) {
	    if($list eq "'") {
		$darktime = sprintf("%8.3f",$junk[1]);
		$list = $arc_file;
		$narcs++;
	    } else {
		$list = join(',',$list,$arc_file);
		$narcs++;
	    }
	}
    }
    if($narcs > 0) {
	$list = $list."'";
	$list =~ s/.fits//g;
	$keyword = 'ARC';
	$value   =  $list;
	$comment = 'arc images';
	$line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
	print LOG $line,"\n";
#
# now the corresponding darks
#
	@darks = ();
	for($j = 0 ; $j <= $#files ; $j++) {
	    if($files[$j] =~ m/dark/) {
		push(@darks,$files[$j]);
	    }
	}
	$list = "'".$darks[0];
	for($j = 1; $j <= $#darks; $j++) {
	    $list = join(',',$list,$darks[$j]);
	}
	$list = $list."'";
	$list =~ s/.fits//g;
	$keyword = 'DARKARC';
	$len  = length($list);
	print "$keyword length is $len\n";
	if($len <= 69) {
	    $value = $list;
	    $line = sprintf("%-8s= %-70s",$keyword, $value);
	    print LOG $line,"\n";
	} else {
# break into two lines
	    $value   =  substr($list,0,68)."&'";
	    $line2   =  "'".substr($list,68);
	    $comment = 'dark frames for arcs';
	    $line = sprintf("%-8s= %-70s",$keyword, $value);
	    print LOG $line,"\n";
#
	    $keyword = 'CONTINUE';
	    $line = sprintf("%-8s %-70s",$keyword, $line2);
	    print LOG $line,"\n";
	    
	    $keyword = 'LONGSTRN';
	    $value   = '';
	    $line = sprintf("%-9s",$keyword);
	    print LOG $line,"\n";
	}
    }
#
#--------------------------------------------------------------------
#
# Reduction steps
#
    $line = 'COMMENT   process the following reduction steps (yes=1/no=0)';
    print LOG $line,"\n";
#
#  ramp fitting and dark subtraction
#
    $keyword = 'S01PROC';
    $value   = 1;
    $comment = 'up-the-ramp fitting (if needed) and dark subtraction';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# distortion map creation
#
    $keyword = 'S02PROC';
    if($science_type[$file_in_list] =~ m/open/) {
	$value = 0;
    } else {
	$value  = 1;
    }
    $comment = 'distortion map creation';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# flat-fielding and 2D slit extraction
#
    $keyword = 'S03PROC';
    if($science_type[$file_in_list] =~ m/open/) {
	$value = 0;
    } else {
	$value  = 1;
    }
    $comment = 'flat-fielding and 2D slit extraction';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# wavelength calibration
#
    $keyword = 'S04PROC';
    if($science_type[$file_in_list] =~ m/open/) {
	$value = 0;
    } else {
	$value  = 1;
    }
    $comment = 'flat-fielding and 2D slit extraction';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# sky-subtraction
#
    $keyword = 'S05PROC';
    if($science_type[$file_in_list] =~ m/open/) {
	$value = 0;
    } else {
	$value  = 1;
    }
    $comment = 'sky-subtraction';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# linearisation
#
    $keyword = 'S06PROC';
    $value   = 1;
    $comment = 'linerisation';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# extraction
#
    $keyword = 'S07PROC';
    if($science_type[$file_in_list] =~ m/open/) {
	$value = 0;
    } else {
	$value  = 1;
    }
    $comment = 'extraction';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# telluric star processing
#
    $keyword = 'S08PROC';
    if($science_type[$file_in_list] =~ m/open/) {
	$value = 0;
    } else {
	$value  = 1;
    }
    $comment = 'telluric star processing';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
# telluric correction
#
    $keyword = 'S09PROC';
    if($science_type[$file_in_list] =~ m/open/) {
	$value = 0;
    } else {
	$value  = 1;
    }
    $comment = 'wavelength calibration';
    $line = sprintf("%-8s= %-40s /%-40s",$keyword, $value, $comment);
    print LOG $line,"\n";
#
    $line = 'END';
    print LOG $line,"\n";
    
    close(LOG);
#    $line = join(',','mmirs_imaging',"'".$logfile."'",'/verbose','/crosstalk','/clean');
    if($two_point == 1) {
	$line = join(',','mmirs_imaging',"'".$logfile."'",'/verbose','/clean','crosstalk='.$crosstalk,'/two_point');
    } else {
	$line = join(',','mmirs_imaging',"'".$logfile."'",'/verbose','/clean','crosstalk='.$crosstalk);
    }
    print IDL $line,"\n";
}
close(IDL);

print "idl_batch is $idl_batch\n";

#
#--------------------------------------------------------------------------
#
#sub print_header_line{
#    my($keyword, $value, $comment) =@_;
    

#--------------------------------------------------------------------------
# check CFITSIO status
#
sub check_status {
    my ($s) = @_;
#    print "check_status $s\n";
    if ($s != 0) {
        my ($txt);
        Astro::FITS::CFITSIO::fits_get_errstatus($s,$txt);
        carp "CFITSIO error: $txt";
        return 0;
    }
    return 1;
}
#
