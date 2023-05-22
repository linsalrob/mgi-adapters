use strict;
use Rob;
my $rob = new Rob;

=pod

This is a simple script that reads the trimming_slurm directory and then summarises how many samples we trimmed

=cut


my $dir = shift;
if (!$dir && -e "trimming_slurm") {
	print STDERR "Using the directory trimming slurm\n";
	$dir = "trimming_slurm";
}
if (!$dir) {
	die "You didn't provide a directory and we can't find the default trimming_slurm directory\n";
}

my ($r1tot, $r2tot, $r1trim, $r2trim, $seqs) = (0, 0, 0, 0, 0);


opendir(DIR, $dir) || die "can't open $dir";
foreach my $file (grep {$_ !~ /^\./} readdir(DIR)) {
	open(IN, "$dir/$file") || die "Can't open $dir/$file";
	while (<IN>) {
		if (m/Total sequences: R1 (\d+) R2 (\d+)/) {
			$r1tot+=$1;
			$r2tot+=$2;
			$seqs++;
		}
		if (m/Sequences trimmed: R1 (\d+) R2 (\d+)/) {
			$r1trim+=$1;
			$r2trim+=$2;
		}
	}
}
print "Directory\tSequence files\tR1 total\tR2 total\tR1 trimmed\tR2 trimmed\n";
print join("\t", $dir, $rob->commify($seqs), $rob->commify($r1tot), $rob->commify($r2tot), $rob->commify($r1trim), $rob->commify($r2trim)), "\n";

