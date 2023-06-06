#!/usr/bin/perl
open GO,"$ARGV[0]";
open GOA,"$ARGV[1]";
open OUT,">$ARGV[2]";

while(<GO>)
{
        chomp;
        
        
        if($_ =~ /\A\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+.*)\Z/)
        {
               
               $a{$3}=$2;
               $d{$3}=$1;
               #print "$2\t";
                
        }

}


while(<GOA>)
{
        chomp;
        
        #s/ /_/g;
        if($_ =~ /\A(.*)\s+(\d+)\Z/)
        {
                $z=$1;
               #print OUT "$_\t";
               @x = split /\|/, $1;

               $id = 1;
               foreach $x (@x){
               
                    if($x=~/\S_(.*)\Z/)
                    {

                       $qq=$1;
                       #print "$1\t";
                      #print OUT "$a{$1}"."|";
                       if($id ==1){
                           $td=$a{$1};
                       }

                       if(($id !=1)&&(exists$a{$1})){
                           $td= $td."\|".$a{$1};
                       }

                       $id ++;
                     }
               }
               
           print OUT "$a{$qq}\t$z\t$td\t$d{$qq}\n";     
        }

}


