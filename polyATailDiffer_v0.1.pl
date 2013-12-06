#!/usr/bin/perl/ -w
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use List::Util qw ( sum );
use Storable;
use Data::Dumper::Names;

######################################################################################################################################################
#
#	Description
#		This is a perl script to find the differences of the polyA tails between multiple samples. It reads a plain txt file specifying the path of the polyA SAM file 
#	and their sample names, as well as part of the output of the polyA tail discoverer.
#
#	Input
#		--polyASAMPathFile=				file path; a plain txt file specifying the path of the polyA SAM file and their sample names; mandatory;
#		--geneBasedPolyAInfoPath=		file path; path of geneBasedEnd3SensePolyA.info.txt output from polyATailDiscoverer;
#		--polyAClstrInfoPath=			file path; path of clusterInfo.txt output from polyATailDiscoverer;
#		--outDir=						directory for output;
#
#	Usage
#		
#		perl polyATailDiffer_v0.1.pl --polyASAMPathFile=./HM1RhmPairEndPolyASAM.txt --polyAClstrInfoPath=/Volumes/B_MPro2TB/NGS/results/EHI_polyA_pairEndOnly_5xA/finalSAM/dynPct70_uNum10.8_dNum10.8_uTrk0.0_dTrk5.5_Len18.999_NM99_NH100_rmRdno/R1P2A8G20C0B200S200K6.6P50.50.0.50/clusterInfo.log.txt --geneBasedPolyAInfoPath=/Volumes/B_MPro2TB/NGS/results/EHI_polyA_pairEndOnly_5xA/finalSAM/dynPct70_uNum10.8_dNum10.8_uTrk0.0_dTrk5.5_Len18.999_NM99_NH100_rmRdno/R1P2A8G20C0B200S200K6.6P50.50.0.50/geneBasedEnd3SensePolyA.info.txt --outDir=./polyATailDiscoverer
#
#	Assumption
#
#	History:
#		
#		v0.1
#		-built based on polyATailDiscover v0.8;
#
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#
#----------Read parameters ----------#
use vars qw ($polyASAMPathFile $geneBasedPolyAInfoPath $polyAClstrInfoPath $outDir);
my ($polyASAMPathFile, $geneBasedPolyAInfoPath, $polyAClstrInfoPath, $outDir) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

my $polyASAMPathFileHsh_ref = readPolyASAMPathFile($polyASAMPathFile);
my $allSamplePolyAPileupHsh_ref = readAllPolyASAM($polyASAMPathFileHsh_ref);
my $geneBasedPolyAInfoHsh_ref = readGeneBasedPolyAInfo($geneBasedPolyAInfoPath);
my $polyAClstrInfoHsh_ref = readPolyAClusterInfo($polyAClstrInfoPath);

my @clstrTypeAry =  ('allEnd3SenseClstrID', 'confidentAltClstr', 'allIntrnlClstrID', 'confidentIntrnalSensePolyAClstr');
my ($clstrBasedSamplePolyACountHsh_ref, $sampleTotalPolyACountHsh_ref, $geneBasedSamplePolyACountHsh_ref, $originalWiggleOutputAryHsh_ref, $scaledWiggleOutputAryHsh_ref) = countPolyAInSample($geneBasedPolyAInfoHsh_ref, $polyAClstrInfoHsh_ref, $allSamplePolyAPileupHsh_ref, \@clstrTypeAry);
generateWiggleForPolyASite($originalWiggleOutputAryHsh_ref, "$outDir/wiggle/original/");
generateWiggleForPolyASite($scaledWiggleOutputAryHsh_ref, "$outDir/wiggle/scaled/");

my $clstrTypeForNormalization = 'allEnd3SenseClstrID';
my $minPrprtnShift = 2;
my $minClstrRdNum = 10;

calculateAndPrintPolyAChanges($clstrBasedSamplePolyACountHsh_ref, $sampleTotalPolyACountHsh_ref, $clstrTypeForNormalization, $geneBasedSamplePolyACountHsh_ref, "$outDir/log", $minPrprtnShift, $minClstrRdNum);

exit;
#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	$outDir = "./polyATailDiscoverer";

	foreach my $param (@ARGV) {
		if ($param =~ m/--polyASAMPathFile=/) {$polyASAMPathFile = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--geneBasedPolyAInfoPath=/) {$geneBasedPolyAInfoPath = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--polyAClstrInfoPath=/) {$polyAClstrInfoPath = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);} 
	}
	
	chop $outDir if ($outDir =~ m/\/$/); #---remove the last slash

	die "Cant read polyASAMPathFile $polyASAMPathFile" if not -s $polyASAMPathFile;
	die "Cant read geneBasedPolyAInfoPath $geneBasedPolyAInfoPath" if not -s $geneBasedPolyAInfoPath;
	die "Cant read polyAClstrInfoPath $polyAClstrInfoPath" if not -s $polyAClstrInfoPath;

	my @paramAry = ($polyASAMPathFile, $geneBasedPolyAInfoPath, $polyAClstrInfoPath, $outDir);

	system "mkdir -p -m 777 $outDir/";
	open (PARAM, ">$outDir/parameters.txt");
	print PARAM Dumper(\@paramAry);
	close PARAM;
	
	return @paramAry;
}
########################################################################## readAndCountPolyASAM
sub readAllPolyASAM {
	
	my %polyASAMPathFileHsh = %{$_[0]};
	
	my %allSamplePolyAPileupHsh;
	foreach my $sampleName (sort {$a cmp $b} keys %polyASAMPathFileHsh) {
		foreach my $rep (sort {$a cmp $b} keys %{$polyASAMPathFileHsh{$sampleName}}) {
			my $polyASAMPath = ${$polyASAMPathFileHsh{$sampleName}}{$rep};
			my $polyAPileupHsh_ref = pileupPolyASAM($polyASAMPath);
			my %polyAPileupHsh = %{$polyAPileupHsh_ref};
			
			foreach my $cntg (keys %polyAPileupHsh) {
				foreach my $polyASite (keys %{$polyAPileupHsh{$cntg}}) {
					foreach my $strnd (keys %{${$polyAPileupHsh{$cntg}}{$polyASite}}) {
						${${${${$allSamplePolyAPileupHsh{$sampleName}}{$rep}}{$cntg}}{$polyASite}}{$strnd} = ${${$polyAPileupHsh{$cntg}}{$polyASite}}{$strnd};
					}
				}
			}
		}
	}
	
	return \%allSamplePolyAPileupHsh;
	
	########################################################################## pileupPolyASAM
	sub pileupPolyASAM {

		my $samPathToRead = $_[0];
	
		my %polyAPileupHsh;
		
		#---read the SAM
		open (INFILE, $samPathToRead);
	
		#---define the start time and counters
		my $intervalStart = time();
		my $lineProc = my $progCount = my $totalTimeSpent = 0;
	
		print "\n\n################# Start Reading SAM #################\n\n";
	
		#---get the total line number and determine the interval size
		my ($fileTotalLineNum,  $intervalSize)= checkFileSizeAndDefineIntervalSize($samPathToRead, 1000);
	
		#---define the SAM Table
		my $SAMFlagTableHsh_ref = defineSAMFlagTable();
		my %SAMFlagTableHsh = %{$SAMFlagTableHsh_ref};
	
		#-----print progress bar
		printProgressScale("Summarizing the SAM file", 50);
	
		#---Start reading samPath
		while (my $theLine = <INFILE>) {
	
			#---report the progress
			$lineProc++; $progCount++;
			($progCount, $intervalStart) = reportProgress($progCount, $lineProc, $intervalSize, $fileTotalLineNum, $intervalStart) if ($progCount >= $intervalSize);
	
			next if ($theLine =~ m/^\@/);#---skip the commment and info lines
			
			chomp $theLine;
			my @theLineSplt = split /\t/, $theLine;
		
			#---store all SAM bits in a Hsh
			my $SAMFlag = $theLineSplt[1];
			my $SAMBitStr = $SAMFlagTableHsh{$SAMFlag};
			my @SAMBitAry = split /\+/, $SAMBitStr;
			my %SAMBitHsh;
			foreach my $SAMBit (@SAMBitAry) {$SAMBitHsh{$SAMBit}++;}
	
			if (not exists $SAMBitHsh{4}) {#---matched alignment or end of the file, supposed to be all alinged, as the unaligned reads are removed 
				
				#---Check the current contig
				my $curntCntg = $theLineSplt[2];
				my $length = length $theLineSplt[9];
				
				#---get the basic info of the read
				my $strnd;
				
				if (not exists $SAMBitHsh{16}) {$strnd = "+";} 
				else {$strnd = "-";}
					
				my $start = $theLineSplt[3];
	
				#---Generate an array contained the aligned position based on the cigar string
				my $cigarStr = $theLineSplt[5];
					
				my $end;
				
				if ($cigarStr =~ m/N/) {#---spliced 
					my $genomicLength = 0;
					$genomicLength += $1 while ($cigarStr =~ /(\d+)M/g);
					$genomicLength += $1 while ($cigarStr =~ /(\d+)N/g);
					$end = $start + $genomicLength - 1;
				} else {#---not spliced
					$end = $start + $length - 1;
				}
				
				my $polyASite;
				
				if ($strnd eq "+") {
					$polyASite = $end;
				} else {
					$polyASite = $start;
				}
				
				#---Pileup the reads in the contig specfic hsh
				${${$polyAPileupHsh{$curntCntg}}{$polyASite}}{$strnd}++;
	
			}#---end of if ($theLineSplt[1] != 4)
		}#---end of while (my $theLine = <INFILE>)
	
		print "\n\n################# End Reading SAM #################\n\n";
	
		return \%polyAPileupHsh;
	}
}
########################################################################## defineSAMFlagTable
sub defineSAMFlagTable {
	
#
#copied from http://bioinformatics.bc.edu/chuanglab/wiki/index.php/SAM_pairwise_flag_translation_table
#
#0x0001 1 the read is paired in sequencing, no matter whether it is mapped in a pair 
#0x0002 2 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment)  
#0x0004 4 the query sequence itself is unmapped 
#0x0008 8 the mate is unmapped  
#0x0010 16 strnd of the query (0 for forward; 1 for reverse strnd) 
#0x0020 32 strnd of the mate  
#0x0040 64 the read is the ﬁrst read in a pair  
#0x0080 128 the read is the second read in a pair 
#0x0100 256 the alignment is not primary (a read having split hits may have multiple primary alignment records) 
	
	my %SAMFlagTableHsh = (

		0 => "0",
		1 => "1",
		2 => "2",
		3 => "1+2",
		4 => "0+4",
		5 => "1+4",
		6 => "0+2+4",
		7 => "1+2+4",
		8 => "0+8",
		9 => "1+8",
		10 => "0+2+8",
		11 => "1+2+8",
		12 => "0+4+8",
		13 => "1+4+8",
		14 => "0+2+4+8",
		15 => "1+2+4+8",
		16 => "0+16",
		17 => "1+16",
		18 => "0+2+16",
		19 => "1+2+16",
		20 => "0+4+16",
		21 => "1+4+16",
		22 => "0+2+4+16",
		23 => "1+2+4+16",
		24 => "0+8+16",
		25 => "1+8+16",
		26 => "0+2+8+16",
		27 => "1+2+8+16",
		28 => "0+4+8+16",
		29 => "1+4+8+16",
		30 => "0+2+4+8+16",
		31 => "1+2+4+8+16",
		32 => "0+32",
		33 => "1+32",
		34 => "0+2+32",
		35 => "1+2+32",
		36 => "0+4+32",
		37 => "1+4+32",
		38 => "0+2+4+32",
		39 => "1+2+4+32",
		40 => "0+8+32",
		41 => "1+8+32",
		42 => "0+2+8+32",
		43 => "1+2+8+32",
		44 => "0+4+8+32",
		45 => "1+4+8+32",
		46 => "0+2+4+8+32",
		47 => "1+2+4+8+32",
		48 => "0+16+32",
		49 => "1+16+32",
		50 => "0+2+16+32",
		51 => "1+2+16+32",
		52 => "0+4+16+32",
		53 => "1+4+16+32",
		54 => "0+2+4+16+32",
		55 => "1+2+4+16+32",
		56 => "0+8+16+32",
		57 => "1+8+16+32",
		58 => "0+2+8+16+32",
		59 => "1+2+8+16+32",
		60 => "0+4+8+16+32",
		61 => "1+4+8+16+32",
		62 => "0+2+4+8+16+32",
		63 => "1+2+4+8+16+32",
		64 => "0+64",
		65 => "1+64",
		66 => "0+2+64",
		67 => "1+2+64",
		68 => "0+4+64",
		69 => "1+4+64",
		70 => "0+2+4+64",
		71 => "1+2+4+64",
		72 => "0+8+64",
		73 => "1+8+64",
		74 => "0+2+8+64",
		75 => "1+2+8+64",
		76 => "0+4+8+64",
		77 => "1+4+8+64",
		78 => "0+2+4+8+64",
		79 => "1+2+4+8+64",
		80 => "0+16+64",
		81 => "1+16+64",
		82 => "0+2+16+64",
		83 => "1+2+16+64",
		84 => "0+4+16+64",
		85 => "1+4+16+64",
		86 => "0+2+4+16+64",
		87 => "1+2+4+16+64",
		88 => "0+8+16+64",
		89 => "1+8+16+64",
		90 => "0+2+8+16+64",
		91 => "1+2+8+16+64",
		92 => "0+4+8+16+64",
		93 => "1+4+8+16+64",
		94 => "0+2+4+8+16+64",
		95 => "1+2+4+8+16+64",
		96 => "0+32+64",
		97 => "1+32+64",
		98 => "0+2+32+64",
		99 => "1+2+32+64",
		100 => "0+4+32+64",
		101 => "1+4+32+64",
		102 => "0+2+4+32+64",
		103 => "1+2+4+32+64",
		104 => "0+8+32+64",
		105 => "1+8+32+64",
		106 => "0+2+8+32+64",
		107 => "1+2+8+32+64",
		108 => "0+4+8+32+64",
		109 => "1+4+8+32+64",
		110 => "0+2+4+8+32+64",
		111 => "1+2+4+8+32+64",
		112 => "0+16+32+64",
		113 => "1+16+32+64",
		114 => "0+2+16+32+64",
		115 => "1+2+16+32+64",
		116 => "0+4+16+32+64",
		117 => "1+4+16+32+64",
		118 => "0+2+4+16+32+64",
		119 => "1+2+4+16+32+64",
		120 => "0+8+16+32+64",
		121 => "1+8+16+32+64",
		122 => "0+2+8+16+32+64",
		123 => "1+2+8+16+32+64",
		124 => "0+4+8+16+32+64",
		125 => "1+4+8+16+32+64",
		126 => "0+2+4+8+16+32+64",
		127 => "1+2+4+8+16+32+64",
		128 => "0+128",
		129 => "1+128",
		130 => "0+2+128",
		131 => "1+2+128",
		132 => "0+4+128",
		133 => "1+4+128",
		134 => "0+2+4+128",
		135 => "1+2+4+128",
		136 => "0+8+128",
		137 => "1+8+128",
		138 => "0+2+8+128",
		139 => "1+2+8+128",
		140 => "0+4+8+128",
		141 => "1+4+8+128",
		142 => "0+2+4+8+128",
		143 => "1+2+4+8+128",
		144 => "0+16+128",
		145 => "1+16+128",
		146 => "0+2+16+128",
		147 => "1+2+16+128",
		148 => "0+4+16+128",
		149 => "1+4+16+128",
		150 => "0+2+4+16+128",
		151 => "1+2+4+16+128",
		152 => "0+8+16+128",
		153 => "1+8+16+128",
		154 => "0+2+8+16+128",
		155 => "1+2+8+16+128",
		156 => "0+4+8+16+128",
		157 => "1+4+8+16+128",
		158 => "0+2+4+8+16+128",
		159 => "1+2+4+8+16+128",
		160 => "0+32+128",
		161 => "1+32+128",
		162 => "0+2+32+128",
		163 => "1+2+32+128",
		164 => "0+4+32+128",
		165 => "1+4+32+128",
		166 => "0+2+4+32+128",
		167 => "1+2+4+32+128",
		168 => "0+8+32+128",
		169 => "1+8+32+128",
		170 => "0+2+8+32+128",
		171 => "1+2+8+32+128",
		172 => "0+4+8+32+128",
		173 => "1+4+8+32+128",
		174 => "0+2+4+8+32+128",
		175 => "1+2+4+8+32+128",
		176 => "0+16+32+128",
		177 => "1+16+32+128",
		178 => "0+2+16+32+128",
		179 => "1+2+16+32+128",
		180 => "0+4+16+32+128",
		181 => "1+4+16+32+128",
		182 => "0+2+4+16+32+128",
		183 => "1+2+4+16+32+128",
		184 => "0+8+16+32+128",
		185 => "1+8+16+32+128",
		186 => "0+2+8+16+32+128",
		187 => "1+2+8+16+32+128",
		188 => "0+4+8+16+32+128",
		189 => "1+4+8+16+32+128",
		190 => "0+2+4+8+16+32+128",
		191 => "1+2+4+8+16+32+128",
		192 => "0+64+128",
		193 => "1+64+128",
		194 => "0+2+64+128",
		195 => "1+2+64+128",
		196 => "0+4+64+128",
		197 => "1+4+64+128",
		198 => "0+2+4+64+128",
		199 => "1+2+4+64+128",
		200 => "0+8+64+128",
		201 => "1+8+64+128",
		202 => "0+2+8+64+128",
		203 => "1+2+8+64+128",
		204 => "0+4+8+64+128",
		205 => "1+4+8+64+128",
		206 => "0+2+4+8+64+128",
		207 => "1+2+4+8+64+128",
		208 => "0+16+64+128",
		209 => "1+16+64+128",
		210 => "0+2+16+64+128",
		211 => "1+2+16+64+128",
		212 => "0+4+16+64+128",
		213 => "1+4+16+64+128",
		214 => "0+2+4+16+64+128",
		215 => "1+2+4+16+64+128",
		216 => "0+8+16+64+128",
		217 => "1+8+16+64+128",
		218 => "0+2+8+16+64+128",
		219 => "1+2+8+16+64+128",
		220 => "0+4+8+16+64+128",
		221 => "1+4+8+16+64+128",
		222 => "0+2+4+8+16+64+128",
		223 => "1+2+4+8+16+64+128",
		224 => "0+32+64+128",
		225 => "1+32+64+128",
		226 => "0+2+32+64+128",
		227 => "1+2+32+64+128",
		228 => "0+4+32+64+128",
		229 => "1+4+32+64+128",
		230 => "0+2+4+32+64+128",
		231 => "1+2+4+32+64+128",
		232 => "0+8+32+64+128",
		233 => "1+8+32+64+128",
		234 => "0+2+8+32+64+128",
		235 => "1+2+8+32+64+128",
		236 => "0+4+8+32+64+128",
		237 => "1+4+8+32+64+128",
		238 => "0+2+4+8+32+64+128",
		239 => "1+2+4+8+32+64+128",
		240 => "0+16+32+64+128",
		241 => "1+16+32+64+128",
		242 => "0+2+16+32+64+128",
		243 => "1+2+16+32+64+128",
		244 => "0+4+16+32+64+128",
		245 => "1+4+16+32+64+128",
		246 => "0+2+4+16+32+64+128",
		247 => "1+2+4+16+32+64+128",
		248 => "0+8+16+32+64+128",
		249 => "1+8+16+32+64+128",
		250 => "0+2+8+16+32+64+128",
		251 => "1+2+8+16+32+64+128",
		252 => "0+4+8+16+32+64+128",
		253 => "1+4+8+16+32+64+128",
		254 => "0+2+4+8+16+32+64+128",
		255 => "1+2+4+8+16+32+64+128",
		256 => "256+0",
		272 => "256+0+16",
	);
	 
	 return \%SAMFlagTableHsh;
}
########################################################################## checkFileSizeAndDefineIntervalSize
sub checkFileSizeAndDefineIntervalSize {
    
    my $fileToCheckPath = $_[0];
    my $linesToSample = $_[1];
	
	#---make sure $linesToSample is a non-zero number, if not set to 10000
	my $linesToSampleInt = int $linesToSample;
	$linesToSampleInt = 100000 if (($linesToSampleInt != $linesToSample) or ($linesToSampleInt == 0));

    #---get the filename from the path
    my @fileToCheckPathSplt = split /\//, $fileToCheckPath;
    my $fileToCheckName = $fileToCheckPathSplt[-1];
    
    print "Estimating the number of lines in $fileToCheckName......";
	
	#---estimate the number of lines in the file
	open (INFILE, $fileToCheckPath) || die "Can't open $fileToCheckPath.\n";
	my $tmpFilePath = $fileToCheckPath."_tmp.txt";
	my $fileToCheckSize = -s "$fileToCheckPath";
	my $fileToCheckSizeTotalLineNum = 0;
	if ($fileToCheckSize < 100000000) {#----small file <1000MB
		chomp ($fileToCheckSizeTotalLineNum = `wc -l < $fileToCheckPath`);
	} else {
		system "tail -$linesToSampleInt $fileToCheckPath >$tmpFilePath";
		my $tmpFileSize = -s "$tmpFilePath";
		system "rm $tmpFilePath";
		my $fileToCheckSizeTotalLineNum = int (($fileToCheckSize/$tmpFileSize)*$linesToSampleInt);
	}
	print ":Estimated ".$fileToCheckSizeTotalLineNum." lines.\n";

	my $intervalSize = int ($fileToCheckSizeTotalLineNum/100); #---define as 
	$intervalSize = 1000000 if  ($intervalSize > 1000000);
	
	return ($fileToCheckSizeTotalLineNum, $intervalSize);

}
########################################################################## reportProgress
sub reportProgress {

	my $progCount = $_[0];
	my $lineProc = $_[1];
	my $intervalSize = $_[2];
	my $fileTotalLineNum = $_[3];
	my $intervalStart = $_[4];

	$progCount=0;
	my $intervalEnd = time();
	my $timeElapsed = $intervalEnd - $intervalStart;
	$timeElapsed = sprintf ("%.2f", $timeElapsed);
	my $estimatedEnd = (($fileTotalLineNum - $lineProc)*$timeElapsed)/$intervalSize;
	$estimatedEnd = sprintf ("%.2f", $estimatedEnd/60);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	updateProgressBar("End in $estimatedEnd min", $lineProc, $fileTotalLineNum, 50, 5);
	$intervalStart = time();
	
	return ($progCount, $intervalStart);
	
	########################################################################## updateProgressBar
	sub updateProgressBar {
	
		my $strToPrint = $_[0];
		my $progressCount = $_[1];
		my $totalCount = $_[2];
		my $scaleMax = $_[3];
		my $extraWhiteSpaceNum = $_[4]; #---erase the longer infos during the progress
	
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	

		my $progressPct = int (($progressCount/$totalCount)*$scaleMax);	

		my $progressBar = "|";
		for my $i (1..$progressPct) {$progressBar .= ">";}
		for my $i (($progressPct+1)..$scaleMax) {$progressBar .= " ";}
		$progressBar .= "|";

		my $extraWhiteSpaceStr = "";
		for my $i (1..$extraWhiteSpaceNum) {$extraWhiteSpaceStr .= " ";}
	
		print $progressBar."[$runTime]".$strToPrint.$extraWhiteSpaceStr."\r";
	}
}
########################################################################## printProgressScale
sub printProgressScale {

	my $strToPrint = $_[0];
	my $scaleMax = $_[1];

	my $scaleSpace = "|";
	for my $i (1..$scaleMax) {$scaleSpace .= "-";}
	$scaleSpace .= "|100%";
	
	print $strToPrint."\n";
	print $scaleSpace."\n";
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
########################################################################## generateWiggleForPolyASite
sub generateWiggleForPolyASite {
	
	#---in/out: generateWiggleForPolyASite($wiggleOutputAryHsh_ref, $subOutDir);
	
	my %wiggleOutputAryHsh = %{$_[0]}; #---push @{${${${$wiggleOutputAryHsh{$clstrType}}{$sampleName}}{'pooled'}}{$strnd}}, $wiggleLine;
	my $subOutDir = $_[1];

	system "mkdir -p -m 777 $subOutDir";
	
	print "Printing wiggles\n";

	my %pooledAllClustrWiggleHsh;

	foreach my $clstrType (keys %wiggleOutputAryHsh) {
		foreach my $sampleName (keys %{$wiggleOutputAryHsh{$clstrType}}) {
			foreach my $rep (keys %{${$wiggleOutputAryHsh{$clstrType}}{$sampleName}}) {
				foreach my $strnd (keys %{${${$wiggleOutputAryHsh{$clstrType}}{$sampleName}}{$rep}}) {
					open (WIGGLE, ">$subOutDir/$clstrType.$sampleName.$rep.$strnd.wig");
					foreach my $cntg (sort {$a cmp $b} keys %{${${${$wiggleOutputAryHsh{$clstrType}}{$sampleName}}{$rep}}{$strnd}}) {
						foreach my $adjustedSite (sort {$a cmp $b} keys %{${${${${$wiggleOutputAryHsh{$clstrType}}{$sampleName}}{$rep}}{$strnd}}{$cntg}}) {
							my $rdNum = ${${${${${$wiggleOutputAryHsh{$clstrType}}{$sampleName}}{$rep}}{$strnd}}{$cntg}}{$adjustedSite};
							${${${${${$pooledAllClustrWiggleHsh{'allClstr'}}{$sampleName}}{$rep}}{$strnd}}{$cntg}}{$adjustedSite} = $rdNum;
							print WIGGLE $cntg."\t".$adjustedSite."\t".$adjustedSite."\t".$rdNum."\n";
						}
					}
					close WIGGLE;
				}
			}
		}
	}
	
	foreach my $clstrType (keys %pooledAllClustrWiggleHsh) {
		foreach my $sampleName (keys %{$pooledAllClustrWiggleHsh{$clstrType}}) {
			foreach my $rep (keys %{${$pooledAllClustrWiggleHsh{$clstrType}}{$sampleName}}) {
				foreach my $strnd (keys %{${${$pooledAllClustrWiggleHsh{$clstrType}}{$sampleName}}{$rep}}) {
					open (WIGGLE, ">$subOutDir/$clstrType.$sampleName.$rep.$strnd.wig");
					foreach my $cntg (sort {$a cmp $b} keys %{${${${$pooledAllClustrWiggleHsh{$clstrType}}{$sampleName}}{$rep}}{$strnd}}) {
						foreach my $adjustedSite (sort {$a cmp $b} keys %{${${${${$pooledAllClustrWiggleHsh{$clstrType}}{$sampleName}}{$rep}}{$strnd}}{$cntg}}) {
							my $rdNum = ${${${${${$pooledAllClustrWiggleHsh{$clstrType}}{$sampleName}}{$rep}}{$strnd}}{$cntg}}{$adjustedSite};
							print WIGGLE $cntg."\t".$adjustedSite."\t".$adjustedSite."\t".$rdNum."\n";
						}
					}
					close WIGGLE;
				}
			}
		}
	}
	
	return;
}
########################################################################## readPolyASAMPathFile
sub readPolyASAMPathFile {

	#--- my $polyASAMPathFileHsh_ref = readPolyASAMPathFile($polyASAMPathFile);

	my $polyASAMPathFileToRead = $_[0];
	
	my %polyASAMPathFileHsh;
	
	open (SAMPATHFILE, "$polyASAMPathFileToRead");
	while (my $theLine = <SAMPATHFILE>) {
		next if length $theLine < 3;
		chomp $theLine;
		my ($sampleName, $rep, $polyASAMPath) = split /\t/, $theLine;
		${$polyASAMPathFileHsh{$sampleName}}{$rep} = $polyASAMPath;
		die "Cannot open $polyASAMPath.\n" if not -s $polyASAMPath;
		print "SAM of $sampleName replicate $rep checked.\n";
		
	}
	close SAMPATHFILE;
	
	return \%polyASAMPathFileHsh;
}
########################################################################## readAltPolyAIsoformInfo
sub readGeneBasedPolyAInfo {

	#--- my %geneBasedPolyAInfoHsh = readGeneBasedPolyAInfo($geneBasedPolyAInfoPath);
	
	my $geneBasedPolyAInfoPathToRead = $_[0];

	my %geneBasedPolyAInfoHsh;

	print "\nReading gene based polyA information.\n";
	
	open (ALTPOLYAINFO, "$geneBasedPolyAInfoPathToRead");
	while (my $theLine = <ALTPOLYAINFO>) {
		next if $theLine =~ m/^#/;
		next if length $theLine < 3;
		chomp $theLine;
		my ($gene, $strnd, $cntg, $end3, $end5, $end3PolyAClstrNum, $closestDisFrom3End, $closestPeakPos, $closestCovGapSize, $closestClstrID, $farthestDisFrom3End, $farthestPeakPos, $farthestCovGapSize, $farthestClstrID, $mostRdDisFrom3End, $mostRdPeakPos, $mostRdCovGapSize, $mostRdClstrID, $withinPolyACltrNum, $avgWithinPolyACltrRdNum, $mostRdPolyARdNum, $closestPolyARdNum, $farthestPolyARdNum, $allIntrnlClstrID, $confidentAltPolyANum, $allEnd3SenseClstrID, $confidentAltClstr, $confidentAltUTR3PolyA, $confidentIntrnlSensePolyA, $confidentIntrnalSensePolyAClstr) = split /\t/, $theLine;
		${$geneBasedPolyAInfoHsh{$gene}}{'strnd'} = $strnd;
		${$geneBasedPolyAInfoHsh{$gene}}{'cntg'} = $cntg;
		${$geneBasedPolyAInfoHsh{$gene}}{'end3'} = $end3;
		${$geneBasedPolyAInfoHsh{$gene}}{'end5'} = $end5;
		${$geneBasedPolyAInfoHsh{$gene}}{'end3PolyAClstrNum'} = $end3PolyAClstrNum;
		${$geneBasedPolyAInfoHsh{$gene}}{'mostRdDisFrom3End'} = $mostRdDisFrom3End;
		${$geneBasedPolyAInfoHsh{$gene}}{'closestDisFrom3End'} = $closestDisFrom3End;
		${$geneBasedPolyAInfoHsh{$gene}}{'farthestDisFrom3End'} = $farthestDisFrom3End;
		${$geneBasedPolyAInfoHsh{$gene}}{'confidentAltPolyANum'} = $confidentAltPolyANum;
		${$geneBasedPolyAInfoHsh{$gene}}{'confidentAltUTR3PolyA'} = $confidentAltUTR3PolyA;
		${$geneBasedPolyAInfoHsh{$gene}}{'confidentIntrnlSensePolyA'} = $confidentIntrnlSensePolyA;

		@{${$geneBasedPolyAInfoHsh{$gene}}{'allEnd3SenseClstrID'}} = split /;/, $allEnd3SenseClstrID;
		@{${$geneBasedPolyAInfoHsh{$gene}}{'allIntrnlClstrID'}} = split /;/, $allIntrnlClstrID;
		@{${$geneBasedPolyAInfoHsh{$gene}}{'confidentAltClstr'}} = split /;/, $confidentAltClstr;
		@{${$geneBasedPolyAInfoHsh{$gene}}{'confidentIntrnalSensePolyAClstr'}} = split /;/, $confidentIntrnalSensePolyAClstr;
	}
	close ALTPOLYAINFO;
	
	my $geneStored = keys %geneBasedPolyAInfoHsh;
	
	print "PolyA information of $geneStored genes stored.\n";
	
	return \%geneBasedPolyAInfoHsh;
}
########################################################################## readPolyAClusterInfo
sub readPolyAClusterInfo {

	#--- my $polyAClstrInfoHsh_ref = readPolyAClusterInfo($polyAClstrInfoPath);
	
	my $polyAClstrInfoPathToRead = $_[0];

	my %polyAClstrInfoHsh;

	print "\nReading polyA cluster information.\n";
	
	open (POLYACLSTRINFO, "$polyAClstrInfoPathToRead");
	while (my $theLine = <POLYACLSTRINFO>) {
		next if $theLine =~ m/^#/;
		next if length $theLine < 3;
		chomp $theLine;
		my ($clstrID, $cntg, $strnd, $peakPos, $startPos, $endPos, $posStr, $validity, $peakUpToDownStrmAvgRatio, $peakPolyAID, $peakRdCount, $peakAvgPolyATailLen, $nonPeakRdCountSum, $nonPeakRdCountAvg, $clstrSize, $avgUpToDownStrmAvgRatio, $avgRdCount, $polyASiteNumInClstr, $ophanOrClstr, $geneScene, $geneID, $senseGeneHitWithinCheckRange) = split /\t/, $theLine;
		${$polyAClstrInfoHsh{$clstrID}}{'cntg'} = $cntg;
		${$polyAClstrInfoHsh{$clstrID}}{'strnd'} = $strnd;
		${$polyAClstrInfoHsh{$clstrID}}{'peakPos'} = $peakPos;
		${$polyAClstrInfoHsh{$clstrID}}{'startPos'} = $startPos;
		${$polyAClstrInfoHsh{$clstrID}}{'endPos'} = $endPos;
		${$polyAClstrInfoHsh{$clstrID}}{'validity'} = $validity;
	}
	close POLYACLSTRINFO;
	
	my $clstrStored = keys %polyAClstrInfoHsh;
	
	print "Information of $clstrStored polyA clusters stored.\n";
	
	return \%polyAClstrInfoHsh;
}
########################################################################## readPolyAClusterInfo
sub countPolyAInSample {

	#---my $geneBasedPolyAInfoHsh_ref = countPolyAInSample($geneBasedPolyAInfoHsh_ref, $polyAClstrInfoHsh_ref, $allSamplePolyAPileupHsh_ref);
	my %geneBasedPolyAInfoHsh = %{$_[0]};
	my %polyAClstrInfoHsh = %{$_[1]};
	my %allSamplePolyAPileupHsh = %{$_[2]};
	my @clstrTypeAry = @{$_[3]};
	
	my %originalWiggleOutputAryHsh;
	my %clstrBasedSamplePolyACountHsh;
	my %geneBasedSamplePolyACountHsh;
	my %sampleTotalPolyACountHsh;
	
	foreach my $clstrType (@clstrTypeAry) {
		foreach my $sampleName (sort {$a cmp $b} keys %allSamplePolyAPileupHsh) {
			foreach my $rep (sort {$a cmp $b} keys %{$allSamplePolyAPileupHsh{$sampleName}}) {
				${${${$sampleTotalPolyACountHsh{$sampleName}}{$rep}}{$clstrType}}{'rdNum'} = 0;
			}
		}
	}
	
	foreach my $gene (sort {$a cmp $b} keys %geneBasedPolyAInfoHsh) {
		foreach my $sampleName (sort {$a cmp $b} keys %allSamplePolyAPileupHsh) {
			foreach my $rep (sort {$a cmp $b} keys %{$allSamplePolyAPileupHsh{$sampleName}}) {
				${${$geneBasedSamplePolyACountHsh{$gene}}{$sampleName}}{$rep} = 0;
			}
		}
	}
	
	foreach my $gene (sort {$a cmp $b} keys %geneBasedPolyAInfoHsh) {
		print "Checking sample polyA of $gene..........                  \r";
		foreach my $clstrType (@clstrTypeAry) {

			#---skip the gene with no polyA cluster of this type;
			next if ${${$geneBasedPolyAInfoHsh{$gene}}{$clstrType}}[0] eq 'null';
		
			foreach my $clstrID (@{${$geneBasedPolyAInfoHsh{$gene}}{$clstrType}}) {
				my $startPos = ${$polyAClstrInfoHsh{$clstrID}}{'startPos'};
				my $endPos = ${$polyAClstrInfoHsh{$clstrID}}{'endPos'};
				my $strnd = ${$polyAClstrInfoHsh{$clstrID}}{'strnd'};
				my $cntg = ${$polyAClstrInfoHsh{$clstrID}}{'cntg'};
				foreach my $sampleName (sort {$a cmp $b} keys %allSamplePolyAPileupHsh) {

					my %tmpPooledRepPileupHsh;
					foreach my $pos ($startPos..$endPos) {
						${$tmpPooledRepPileupHsh{$cntg}}{$pos} = 0;
					}
					
					foreach my $rep (sort {$a cmp $b} keys %{$allSamplePolyAPileupHsh{$sampleName}}) {
						my $rdNumInClstr = 0;
						foreach my $pos ($startPos..$endPos) {
							if (exists ${${${${$allSamplePolyAPileupHsh{$sampleName}}{$rep}}{$cntg}}{$pos}}{$strnd}) {
								my $rdNum = ${${${${$allSamplePolyAPileupHsh{$sampleName}}{$rep}}{$cntg}}{$pos}}{$strnd};
								$rdNumInClstr += $rdNum;
								
								#----storing info for wiggle replicates
								${$tmpPooledRepPileupHsh{$cntg}}{$pos} += $rdNum;
								my $adjustedSite = $pos - 1;
								${${${${${$originalWiggleOutputAryHsh{$clstrType}}{$sampleName}}{$rep}}{$strnd}}{$cntg}}{$adjustedSite} = $rdNum;
							}
						}
						${${${$sampleTotalPolyACountHsh{$sampleName}}{$rep}}{$clstrType}}{'rdNum'} += $rdNumInClstr;
						${${${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}{$sampleName}}{$rep} = $rdNumInClstr;
						${${${$geneBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$sampleName}}{$rep} += $rdNumInClstr;
					}

					#----storing info for wiggle sample
					foreach my $pos ($startPos..$endPos) {
						if (${$tmpPooledRepPileupHsh{$cntg}}{$pos} > 0) {
							my $rdNum = ${$tmpPooledRepPileupHsh{$cntg}}{$pos};
							my $adjustedSite = $pos - 1;
							${${${${${$originalWiggleOutputAryHsh{$clstrType}}{$sampleName}}{'pooled'}}{$strnd}}{$cntg}}{$adjustedSite} = $rdNum;
						}
					}
				}
				
			}
		}
	}
	
	print "\n";

	my %scaledWiggleOutputAryHsh;
	
	print "\n";

	#---calculate the factor ratio;
	foreach my $clstrType (@clstrTypeAry) {
		my @totalRepRdNumAry;
		foreach my $sampleName (sort {$a cmp $b} keys %allSamplePolyAPileupHsh) {
			foreach my $rep (sort {$a cmp $b} keys %{$allSamplePolyAPileupHsh{$sampleName}}) {
				push @totalRepRdNumAry, ${${${$sampleTotalPolyACountHsh{$sampleName}}{$rep}}{$clstrType}}{'rdNum'};
			}
		}
		my $sampleNum = keys %allSamplePolyAPileupHsh;
		my $avgRepRdNum = sum(@totalRepRdNumAry)/@totalRepRdNumAry;
		my $avgSampleRdNum = sum(@totalRepRdNumAry)/$sampleNum;
		
		foreach my $sampleName (sort {$a cmp $b} keys %allSamplePolyAPileupHsh) {
			my $sampleRdNum = 0;
			foreach my $rep (sort {$a cmp $b} keys %{$allSamplePolyAPileupHsh{$sampleName}}) {
				my $repRdNum = ${${${$sampleTotalPolyACountHsh{$sampleName}}{$rep}}{$clstrType}}{'rdNum'};
				$sampleRdNum += $repRdNum;
				my $repScaleFactor = $avgRepRdNum/$repRdNum;
				${${${$sampleTotalPolyACountHsh{$sampleName}}{$rep}}{$clstrType}}{'scaleFactor'} = $repScaleFactor;
				print $clstrType."\t".$sampleName."\t".$rep."\t".$repRdNum."\t".$repScaleFactor."\n";
				
				#----scale the rdNum is the originalWiggle
				foreach my $strnd (keys %{${${$originalWiggleOutputAryHsh{$clstrType}}{$sampleName}}{$rep}}) {
					foreach my $cntg (keys %{${${${$originalWiggleOutputAryHsh{$clstrType}}{$sampleName}}{$rep}}{$strnd}}) {
						foreach my $adjustedSite (keys %{${${${${$originalWiggleOutputAryHsh{$clstrType}}{$sampleName}}{$rep}}{$strnd}}{$cntg}}) {
							${${${${${$scaledWiggleOutputAryHsh{$clstrType}}{$sampleName}}{$rep}}{$strnd}}{$cntg}}{$adjustedSite} = $repScaleFactor*${${${${${$originalWiggleOutputAryHsh{$clstrType}}{$sampleName}}{$rep}}{$strnd}}{$cntg}}{$adjustedSite};
						}
					}
				}
			}
			
			#----scale the rdNum is the originalWiggle
			my $sampleScaleFactor = $avgSampleRdNum/$sampleRdNum;
			print $clstrType."\t".$sampleName."\t".'pooled'."\t".$sampleRdNum."\t".$sampleScaleFactor."\n";
			foreach my $strnd (keys %{${${$originalWiggleOutputAryHsh{$clstrType}}{$sampleName}}{'pooled'}}) {
				foreach my $cntg (keys %{${${${$originalWiggleOutputAryHsh{$clstrType}}{$sampleName}}{'pooled'}}{$strnd}}) {
					foreach my $adjustedSite (keys %{${${${${$originalWiggleOutputAryHsh{$clstrType}}{$sampleName}}{'pooled'}}{$strnd}}{$cntg}}) {
						${${${${${$scaledWiggleOutputAryHsh{$clstrType}}{$sampleName}}{'pooled'}}{$strnd}}{$cntg}}{$adjustedSite} = $sampleScaleFactor*${${${${${$originalWiggleOutputAryHsh{$clstrType}}{$sampleName}}{'pooled'}}{$strnd}}{$cntg}}{$adjustedSite};
					}
				}
			}
		}
	}
	
	print "\n";
	
	return \%clstrBasedSamplePolyACountHsh, \%sampleTotalPolyACountHsh,\%geneBasedSamplePolyACountHsh, \%originalWiggleOutputAryHsh, \%scaledWiggleOutputAryHsh;
}
########################################################################## readPolyAClusterInfo
sub calculateAndPrintPolyAChanges {

	#--calculateAndPrintPolyAChanges(\%geneBasedPolyAInfoHsh, \%sampleTotalPolyACountHsh, $clstrTypeForNormalization);

	my %clstrBasedSamplePolyACountHsh = %{$_[0]};
	my %sampleTotalPolyACountHsh = %{$_[1]};
	my $clstrTypeForNormalization = $_[2];
	my %geneBasedSamplePolyACountHsh = %{$_[3]};
	my $subOutDir = $_[4];
	my $minPrprtnShift = $_[5];
	my $minClstrRdNum = $_[6];
	
	my $log2minPrprtnShift = log($minPrprtnShift)/log(2);
	
	system "mkdir -p -m 777 $subOutDir";
	print "Calculating changes in isoform proportion\n";

	my %interSampleRatioChangeHsh;
	my %sampleTotalRdNumHsh;
	
	#---get the total number of read in each gene in each sample
	#---get the clstr based rd num
	my %tmpTotalPolyARdNumHsh;
	my %totalPolyARdNumPerGeneHsh;
	foreach my $clstrType (sort {$a cmp $b} keys %clstrBasedSamplePolyACountHsh) {
		foreach my $gene (sort {$a cmp $b} keys %{$clstrBasedSamplePolyACountHsh{$clstrType}}) {
			foreach my $clstrID (sort {$a cmp $b} keys %{${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}) {
				foreach my $sampleName (sort {$a cmp $b} keys %{${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}) {
					foreach my $rep (sort {$a cmp $b} keys %{${${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}{$sampleName}}) {
						${${${$tmpTotalPolyARdNumHsh{$sampleName}}{$rep}}{$gene}}{$clstrID} = ${${${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}{$sampleName}}{$rep};
						${${$totalPolyARdNumPerGeneHsh{$sampleName}}{$rep}}{$gene} = 0;
					}
				}
			}
		}
	}
	
	#-----sum the clstr based rd num
	foreach my $sampleName (keys %tmpTotalPolyARdNumHsh) {
		foreach my $rep (keys %{$tmpTotalPolyARdNumHsh{$sampleName}}) {
			foreach my $gene (keys %{${$tmpTotalPolyARdNumHsh{$sampleName}}{$rep}}) {
				foreach my $clstrID (keys %{${${$tmpTotalPolyARdNumHsh{$sampleName}}{$rep}}{$gene}}) {
					my $rdNum = ${${${$tmpTotalPolyARdNumHsh{$sampleName}}{$rep}}{$gene}}{$clstrID};
					${${$totalPolyARdNumPerGeneHsh{$sampleName}}{$rep}}{$gene} += $rdNum;
				}
			}
		}
	}
	
	
	foreach my $clstrType (sort {$a cmp $b} keys %clstrBasedSamplePolyACountHsh) {
		open (LOG, ">$subOutDir/$clstrType.polyAFoldChange.log.xls");
		
		#----print the header
		foreach my $gene (sort {$a cmp $b} keys %{$clstrBasedSamplePolyACountHsh{$clstrType}}) {
			foreach my $clstrID (sort {$a cmp $b} keys %{${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}) {
				my @outputAry;
				push @outputAry, 'gene';
				push @outputAry, 'clstrID';
				foreach my $sampleName (sort {$a cmp $b} keys %{${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}) {
					foreach my $rep (sort {$a cmp $b} keys %{${${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}{$sampleName}}) {
						my $scaleFactor = sprintf "%.2f", ${${${$sampleTotalPolyACountHsh{$sampleName}}{$rep}}{$clstrTypeForNormalization}}{'scaleFactor'};
						push @outputAry, "$sampleName.$rep [$scaleFactor]";
					}
					push @outputAry, $sampleName."_ScaledRd";
				}
				foreach my $sampleName (sort {$a cmp $b} keys %{${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}) {
					foreach my $rep (sort {$a cmp $b} keys %{${$geneBasedSamplePolyACountHsh{$gene}}{$sampleName}}) {
						push @outputAry, $sampleName."_".$rep.'_prprtn';
					}
					push @outputAry, $sampleName."_RdNumClstr";
					push @outputAry, $sampleName."_RdNumGene";
					push @outputAry, $sampleName."_prprtn";
				}
				foreach my $refSampleName (sort {$a cmp $b} keys %{${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}) {
					foreach my $qrySampleName (sort {$a cmp $b} keys %{${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}) {
						next if $refSampleName eq $qrySampleName;
						push @outputAry, $refSampleName.".vs.".$qrySampleName.'.prprtnLog2FC';
					}
				}

				push @outputAry, 'totalRdNumInClstr';
				push @outputAry, 'passMinRd';
				push @outputAry, 'prtnShift';
				push @outputAry, 'sampleAtLeast2Reads';

				print LOG join "", ((join "\t", @outputAry), "\n");
				last;
			}
			last;
		}
		
		foreach my $gene (sort {$a cmp $b} keys %{$clstrBasedSamplePolyACountHsh{$clstrType}}) {
			foreach my $clstrID (sort {$a cmp $b} keys %{${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}) {
				my @outputAry;
				push @outputAry, $gene;
				push @outputAry, $clstrID;
				my $sampleAtLeast2Reads = 0;
				foreach my $sampleName (sort {$a cmp $b} keys %{${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}) {
					my $sampleTotalScaleRdNum = 0;
					foreach my $rep (sort {$a cmp $b} keys %{${${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}{$sampleName}}) {
						my $rdNumInClstr = ${${${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}{$sampleName}}{$rep};
						my $scaleFactor = ${${${$sampleTotalPolyACountHsh{$sampleName}}{$rep}}{$clstrTypeForNormalization}}{'scaleFactor'};
						my $scaledRdNumInClstr = sprintf "%.2f", $rdNumInClstr*$scaleFactor;
						push @outputAry, $scaledRdNumInClstr;
						$sampleTotalScaleRdNum += $scaledRdNumInClstr;
						$sampleAtLeast2Reads++ if $rdNumInClstr >= 2;
					}
					push @outputAry, $sampleTotalScaleRdNum;
				}

				my $totalRdNumInClstr = 0;
				foreach my $sampleName (sort {$a cmp $b} keys %{${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}) {
					my $sampleTotalRdNumInGene = 0;
					my $sampleTotalRdNumInClstr = 0;

					foreach my $rep (sort {$a cmp $b} keys %{${$geneBasedSamplePolyACountHsh{$gene}}{$sampleName}}) {
						my $repRdNumInGene = ${${$totalPolyARdNumPerGeneHsh{$sampleName}}{$rep}}{$gene};
						my $rdNumInClstr = ${${${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}{$sampleName}}{$rep};
						my $repPolyAProportionInGene = 'null';
						$repPolyAProportionInGene = sprintf "%.2f", $rdNumInClstr/$repRdNumInGene if $repRdNumInGene > 0;
						$sampleTotalRdNumInGene += $repRdNumInGene;
						$sampleTotalRdNumInClstr += $rdNumInClstr;
						$totalRdNumInClstr += $rdNumInClstr;
						push @outputAry, $repPolyAProportionInGene;
					}
					
					my $samplePolyAProportionInGene = 'null';
					$samplePolyAProportionInGene = sprintf "%.2f", $sampleTotalRdNumInClstr/$sampleTotalRdNumInGene if $sampleTotalRdNumInGene > 0;
					${${${$interSampleRatioChangeHsh{$clstrType}}{$gene}}{$clstrID}}{$sampleName} = $samplePolyAProportionInGene;
					push @outputAry, $sampleTotalRdNumInClstr;
					push @outputAry, $sampleTotalRdNumInGene;
					push @outputAry, $samplePolyAProportionInGene;
				}
				
				my $prprtnShift = 'no';
				foreach my $refSampleName (sort {$a cmp $b} keys %{${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}) {
					foreach my $qrySampleName (sort {$a cmp $b} keys %{${${$clstrBasedSamplePolyACountHsh{$clstrType}}{$gene}}{$clstrID}}) {
						next if $refSampleName eq $qrySampleName;
						my $refPrprtn = ${${${$interSampleRatioChangeHsh{$clstrType}}{$gene}}{$clstrID}}{$refSampleName};
						my $qryPrprtn = ${${${$interSampleRatioChangeHsh{$clstrType}}{$gene}}{$clstrID}}{$qrySampleName};
						
						#---null mean no read on all clustr of the gene
						my $prprtnLog2FC = 'null';
						if (($qryPrprtn ne 'null') and ($refPrprtn ne 'null')) {
							
							if (($qryPrprtn > 0) and ($refPrprtn > 0)) {
								$prprtnLog2FC = log($qryPrprtn/$refPrprtn)/log(2);
							} elsif (($qryPrprtn == 0) and ($refPrprtn > 0)) {
								$prprtnLog2FC = -9999;
							} elsif (($qryPrprtn > 0) and ($refPrprtn == 0)) {
								$prprtnLog2FC = 9999;
							} elsif (($qryPrprtn == 0) and ($refPrprtn == 0)) {
								$prprtnLog2FC = 0;
							}
							
							$prprtnShift = 'yes' if abs($prprtnLog2FC) >= $log2minPrprtnShift;
						}
						push @outputAry, $prprtnLog2FC;
					}
				}

				my $passMinRd = 'no';
				$passMinRd = 'yes' if $totalRdNumInClstr >= $minClstrRdNum;
				push @outputAry, $totalRdNumInClstr;
				push @outputAry, $passMinRd;
				push @outputAry, $prprtnShift;
				push @outputAry, $sampleAtLeast2Reads;

				print LOG join "", ((join "\t", @outputAry), "\n");
			}
		}
		close LOG;
	}
}
