#! /usr/bin/perl

use strict;
use warnings;

#expected usage: perl script.pl <conf>.gro <topol>.top <posres>.itp
#should only be called from "Build-CytC script"
#created by R. G. McAllister March 17, 2015.

#variables to keep file names straight.
my $gro = $ARGV[0];
my $top = $ARGV[1];
my $posres = $ARGV[2];

#variables to store indeces.
my $index1;
my $index2;
my @cys_index;
my @h_index;

#arrays to slurp files.
my @groA;
my @topA;
my @posresA;

#Find hydrogens to remove, and store indeces.
open (my $grofile, "<", $gro);
while (<$grofile>) {
	chomp;
	if ($_ =~ m/^\s*(\d+)(.{5})(.{5})\s*(\d+)\s+(-*[\d.]+)\s+(-*[\d.]+)\s+(-*[\d.]+).*$/) {
		my $resnum = $1;
		my $res = $2;
		my $atom = $3;
		my $index = $4;
		my $xpos = $5;
		my $ypos = $6;
		my $zpos = $7;
		map {$_ =~ y/ //d} $resnum,$res,$atom,$index,$xpos,$ypos,$zpos; #remove spaces that shouldn't be there
		push @groA, [$resnum,$res,$atom,$index,[$xpos,$ypos,$zpos],$_]; #last element is original line, so we can easily pop it back off if unchanged.}
		push @cys_index, $index if ($res eq "CYS" and $atom eq "SG");
		push @h_index, $index if ($res eq "HMCR" and $atom =~ m/HA[BC][12]/ );
	} else {
		push @groA, [$_];
	}
}
close $grofile;

{
my @dists;
for my $int (@h_index) { #compute distances
	for my $sulf (@cys_index) {
		my $sqdist = (($groA[$int+1][4][0] -  $groA[$sulf+1][4][0])**2 + ($groA[$int+1][4][1] -  $groA[$sulf+1][4][1])**2 + ($groA[$int+1][4][2] -  $groA[$sulf+1][4][2])**2);
		push @dists, [$int,$sulf,$sqdist];
	}
}
($index1,$index2) = sort {$a <=> $b} (map {$_->[0]} sort {$a->[2] <=> $b->[2]} @dists)[0..1];
}

#now index1 stores lower index of S-overlapping H, and index2 stores (higher) index of other S-overlapping H.

#We need to modify the atom number in the gro file to eliminate two hydrogens.
$groA[1][0] -= 2;

#next we need to modify index numbers and re-generate gro lines for the file
for my $line (@groA[2..$index1]) { #indeces less than index one are not changed, so just re-generate gro lines
	$line->[5] = sprintf("%5d%-5s%5s%5d%8.3f%8.3f%8.3f",@{$line}[0..3],@{$line->[4]}[0..2]);
}
for my $line (@groA[$index1+2..$index2]) { #indeces between h-atoms must be dercremented by one.
	if ($line->[3] == $index1+1) {
		substr($line->[2],-1,1) = "";
	}
	--($line->[3]);
	$line->[5] = sprintf("%5d%-5s%5s%5d%8.3f%8.3f%8.3f",@{$line}[0..3],@{$line->[4]}[0..2]);
}
for my $line (@groA[$index2+2..$groA[1][0]+3]) { #indeces of atoms after second index must be decremented by 2.
	if ($line->[3] == $index2+1) {
		substr($line->[2],-1,1) = "";
	}
	$line->[3] -= 2;
	$line->[5] = sprintf("%5d%-5s%5s%5d%8.3f%8.3f%8.3f",@{$line}[0..3],@{$line->[4]}[0..2]);
}

for my $line (@groA) {
	$line = pop @{$line};
}
splice @groA, $index2+1, 1; #remove the second line first, so that first index is unaffected.
splice @groA, $index1+1, 1;
push @groA,"";

#print the new gro file over top of the old one (but keep the old one for now, for testing purposes)
open(my $groOut, ">", $gro);
print {$groOut} join("\n",@groA);
close $groOut;

#we're done with the gro file. now, so remove the array.
undef @groA;

#now we will work with posres file (because it's easier than the top file.
open(my $posfile, "<", $posres);
while (<$posfile>) {
	push @posresA, $_;
}
close $posfile;

for my $line (@posresA) {
	next unless ($line =~ m/^\s*(\d+)\s*\d+\s*\d+\s*\d+\s*\d+/);
	my $int = $1;
	next if ($int <= $index1);
	--$int;
	--$int if ($int >= $index2);
	$int = sprintf("%6d",$int);
	$line =~ s/^.{6}/$int/;
}

open(my $posOut, ">", $posres); #write over old file (except for testing purposes)
print {$posOut} join('',@posresA);
close $posOut;

#clear out the posres array
undef @posresA;

#finally, we will deal with the topology file
open(my $topfile, "<", "$top");
while (<$topfile>) {
	push @topA, $_;
}
close $topfile;

my @improper = (0,0,0,0); #indeces for hshc-hmcr improper NE2-CD2-CE1-FE
my $imp_done = 0;
my $hshc = 0;
my $hshc_done = 0;
my $hem = 0;
my $hem_done = 0;
my $section = "header"; #easiest to process each section differently
for my $line (@topA) {
	if ($line =~ m/^\[\s+(\w+)/) { #identify what section we're in
		$section = $1;
		next;
	}
	next if ($line =~ m/^#/ or $line =~ m/^$/); #skip directives and blank lines
	next if ($line =~ m/^;/ and not ($line =~ m/HIS rtp HSHC/ or $line =~ m/HMCR rtp HMCR/)); #skip comments except for HSHC and HMCR lines
	next if ($section eq "moleculetype" or $section eq "system" or $section eq "molecules" or $section eq "position_restraints"); #nothing to do in these sections
	if ($section eq "atoms") {
		$line =~ m/\s+(\d+)\s+(\w+)\s+(\d+)\s+(\w+)\s+(\w+)\s+(\d+)(.+\n)/;
		my $index_num = $1;
		my $ff_atom = $2;
		my $res_num = $3;
		my $res_name = $4;
		my $gro_atom = $5;
		my $cg_num = $6;
		my $rest = $7;
		if ($line =~ m/HIS rtp HSHC/) {
			$hshc = 1;
			next;
		}
		if ($hshc and not $hshc_done) {
			$improper[0] = $index_num if ($gro_atom eq "NE2");
			$improper[1] = $index_num if ($gro_atom eq "CD2");
			$improper[2] = $index_num if ($gro_atom eq "CE1");
			++$hshc_done if ((grep {$_} @improper[0..2]) == 3);
		}
		if ($line =~ m/HMCR rtp HMCR/) {
			$hem = 1;
			next;
		}
		if ($hem and not $hem_done) {
			$improper[3] = $index_num if ($gro_atom eq "FE");
			++$hem_done if ($improper[3]);
		}
		next if ($index_num < $index1); #we don't need to modify these
		$line = "" if ($index_num == $index1); #we want to delete this line
		if ($index_num > $index1 and $index_num < $index2) { #decrement index and cgnr. Also fix name of seond H.
			substr($gro_atom,-1,1) = "" if ($index_num == $index1+1);
			--$index_num;
			--$cg_num;
			$line = sprintf("%6d%11s%7d%7s%7s%7d",$index_num,$ff_atom,$res_num,$res_name,$gro_atom,$cg_num) . $rest;
			next;
		}
		$line = "" if ($index_num == $index2); #we want to delete this line
		if ($index_num > $index2) { #double-decrement index and cgnr. Also fix name of seond H.
			substr($gro_atom,-1,1) = "" if ($index_num == $index2+1);
			$index_num -= 2;
			$cg_num -= 2;
			$line = sprintf("%6d%11s%7d%7s%7s%7d",$index_num,$ff_atom,$res_num,$res_name,$gro_atom,$cg_num) . $rest;
			next;
		}
	}
	if ($section eq "bonds" or $section eq "pairs") {
		$line =~ m/\s*(\d+)\s+(\d+)(.*\n)/;
		my @nums = ($1,$2);
		my $rest = $3;
		next if ($nums[0] < $index1 and $nums[1] < $index1); #no changes to these lines
		if (grep {$_ == $index1} @nums or grep {$_ == $index2} @nums) { #eliminate these lines
			$line = "";
			next;
		}
		@nums = map {$_ > $index2 ? --$_ : $_} @nums; #reduce values if necessary
		@nums = map {$_ > $index1 ? --$_ : $_} @nums;
		$line = sprintf("%5d%6d",@nums) . $rest;
	}
	if ($section eq "angles") {
		$line =~ m/\s*(\d+)\s+(\d+)\s+(\d+)(.*\n)/;
		my @nums = ($1,$2,$3);
		my $rest = $4;
		next if (!(grep {$_ >= $index1} @nums)); #no changes to these lines
		if (grep {$_ == $index1} @nums or grep {$_ == $index2} @nums) { #eliminate these lines
			$line = "";
			next;
		}
		@nums = map {$_ > $index2 ? --$_ : $_} @nums; #reduce values if necessary
		@nums = map {$_ > $index1 ? --$_ : $_} @nums;
		$line = sprintf("%5d%6d%6d",@nums) . $rest;
	}
	if ($section eq "dihedrals") {
		$line =~ m/\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)(.*\n)/;
		my @nums = ($1,$2,$3,$4);
		my $rest = $5;
		if ($rest =~ m/2/ and not $imp_done) {#first line of the impropers section
			@improper = map {$_ > $index2 ? --$_ : $_} @improper; #reduce values if necessary
			@improper = map {$_ > $index1 ? --$_ : $_} @improper;
			if (!(grep {$_ >= $index1} @nums)) {
				$line = sprintf("%5d%6d%6d%6d",@improper) . "     2 " . "\n" . $line;
				++$imp_done;
			}
		}
		next if (!(grep {$_ >= $index1} @nums)); #no changes to these lines
		if (grep {$_ == $index1} @nums or grep {$_ == $index2} @nums) { #eliminate these lines
			$line = "";
			next;
		}
		@nums = map {$_ > $index2 ? --$_ : $_} @nums; #reduce values if necessary
		@nums = map {$_ > $index1 ? --$_ : $_} @nums;
		$line = sprintf("%5d%6d%6d%6d",@nums) . $rest;
	}
	if ($section eq "cmap") {
		$line =~ m/\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)(.*\n)/;
		my @nums = ($1,$2,$3,$4,$5);
		my $rest = $6;
		next if (!(grep {$_ >= $index1} @nums)); #no changes to these lines
		if (grep {$_ == $index1} @nums or grep {$_ == $index2} @nums) { #eliminate these lines
			$line = "";
			next;
		}
		@nums = map {$_ > $index2 ? --$_ : $_} @nums; #reduce values if necessary
		@nums = map {$_ > $index1 ? --$_ : $_} @nums;
		$line = sprintf("%5d%6d%6d%6d%6d",@nums) . $rest;
	}
}

#print topology over initial file (except for testing purposes)
open(my $topOut, ">", $top);
print {$topOut} join('',@topA);
close $topOut;