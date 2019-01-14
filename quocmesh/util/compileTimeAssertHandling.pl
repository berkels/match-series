#!/usr/bin/perl -w

while($line = <>) {
    if($line =~ /invalid use of incomplete type .*struct AssertAtCompileTime.*/) {
        chomp $line;

        $line .= " This means that a compile time assertion (also known as static assertion) has failed!\n";
    }

    print $line;
}
