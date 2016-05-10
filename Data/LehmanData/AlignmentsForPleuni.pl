
open (IN, "Reference_RT_HIV1.txt");
while (<IN>) {
    />/ && next;
    /\S+/;
    $ref .= $&;
}




foreach $file (@ARGV) {
    #  NEED TO DO ALIGN!
    system "bwa bwasw Reference_RT_HIV1.txt $file > $file.bwa";
    
    open (BWA, "$file.bwa");
    open (OUT, ">$file.aligned");
    $junk = <BWA>;
    while (<BWA>) {
	($n,$pos,$align,$seq) = (split)[0,3,5,9];
	if ($align =~ s/^(\d+)S//) {
	    $s = $1;
	    $seq =~ s/^.{$s}//;
	    $pos += $s;
	}
	if ($align =~ s/(\d+)S$//) {
	    $s = $1;
	    $seq =~ s/.{$s}$//;
	}
	
	open (IN, ">temp");
	print IN ">1\n$ref\n>2\n$seq\n";
	close IN;
	system "clustalw2 temp -output=gde";
	open (IN, "temp.gde");
	{ local $/ = "#";
	  $j = <IN>;
	  $R = <IN>;
	  $R =~ s/[^a-z\-]+//g;
	  $S = <IN>;
	  $S =~ s/[^a-z\-]+//g;
	  $l = length $S;
	  $t = "$R $S ";
	  $t =~ tr/a-z/A-Z/;

	  $tt = join ("", $t =~ /(?=(.).{$l}(.).*( ))/g);
	  $tt =~ s/ (?=\-)//g;

	  $ttt = join ("", $tt =~ /(\S+) (\S+) (\S+ )/g);
	  @cod = $ttt =~ /\S+/g;
	  for $i (0..$#cod) {
	      $gap[$i] = $cod[$i] =~ /\-/;
	  }
	  if ($gap[0] || $gap[1]) {$cod[0] = "NNNNNN"};
	  for $i (1..$#cod) {
	      if ($gap[$i-1] || $gap[$i] || $gap[$i+1]) {
		  $cod[$i] = "NNNNNN";
	      }
	  }
	  $ttt = join (" ", @cod, "");

	  
	  $al = join "", $ttt =~ /\S(\S)\S(\S)\S(\S )/g;

	  @c = $al =~ /\S+/g;
	  for ($i = 0; $c[$i] =~ /(NNN|---)/; $i++) {
	      $c[$i] = "---";
	  }
	  for ($i = $#c; $c[$i] =~ /(NNN|---)/; $i--) {
	      $c[$i] = "---";
	  }
	
  
	  print OUT ">$n\n";
	  print OUT @c;
	  print OUT "\n";
	  
	}


      
    }

	  
}

