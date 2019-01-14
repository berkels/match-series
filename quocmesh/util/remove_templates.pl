#!/usr/bin/perl -w

# truncates template parameters to one char in warning and error messages

while ( defined ($str=<>)) {
  while ( $str =~ /<[^<^>]+>/ ) {
    $str =~ s/<([^<^>])[^<^>]*>/\@$1\@/g;
  }
  $str =~ s/\@(.)\@/<$1>/g;
  $str =~ s/aol:://g;
  $str =~ s/std:://g;
  print $str;
}
