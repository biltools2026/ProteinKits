#!/usr/bin/perl -w
use strict;
my($window,$line,@fileList,$fnum,@label,@ss1,@solve,@wds,$i,$k,$len,$position,$m,$n);
my(@trains_lengths,@tests_lengths,@train_fileList,@test_fileList,$j,$count);
my(@layer1_predict_train,@layer1_predict_test);
my($len1,$len2,$len3,$len4);
my(@aa,%spx_hash);

 %spx_hash=(
        "R"=>"0.105   0.373   0.466  -0.900   0.900   0.528  -0.371" ,
        "K"=>"-0.088   0.066   0.163  -0.889   0.727   0.279  -0.265" ,
        "D"=>"-0.213  -0.417  -0.281  -0.767  -0.900  -0.155  -0.635" ,
        "E"=>"-0.230  -0.241  -0.058  -0.696  -0.868   0.900  -0.582" ,
        "N"=>"-0.213  -0.329  -0.243  -0.674  -0.075  -0.403  -0.529" ,
        "Q"=>"-0.230  -0.110  -0.020  -0.464  -0.276   0.528  -0.371" ,
        "H"=>"0.384   0.110   0.138  -0.271   0.195  -0.031  -0.106" ,
        "Y"=>"0.363   0.417   0.541   0.188  -0.274  -0.155   0.476" ,
        "W"=>"0.479   0.900   0.900   0.900  -0.209   0.279   0.529" ,
        "S"=>"-0.337  -0.637  -0.544  -0.364  -0.265  -0.466  -0.212" ,
        "T"=>"0.402  -0.417  -0.321  -0.199  -0.288  -0.403   0.212" ,
        "G"=>"-0.900  -0.900  -0.900  -0.342  -0.179  -0.900  -0.900" ,
        "P"=>"0.247  -0.900  -0.294   0.055  -0.010  -0.900   0.106" ,
        "A"=>"-0.350  -0.680  -0.677  -0.171  -0.170   0.900  -0.476" ,
        "M"=>"0.110   0.066   0.087   0.337  -0.262   0.652  -0.001" ,
        "C"=>"-0.140  -0.329  -0.359   0.508  -0.114  -0.652   0.476" ,
        "F"=>"0.363   0.373   0.412   0.646  -0.272   0.155   0.318" ,
        "L"=>"0.213  -0.066  -0.009   0.596  -0.186   0.714  -0.053" ,
        "V"=>"0.677  -0.285  -0.232   0.331  -0.191  -0.031   0.900" ,
        "I"=>"0.900  -0.066  -0.009   0.652  -0.186   0.155   0.688" ,
);
$window=5;
@fileList=();
open(IN,"training_dataset.txt")||die "can not open";
    while($line=<IN>){
    	chomp($line);
    	push(@fileList,$line);
    }
close IN || die "can not close";
$count = 0;
open(OUT,">train_svm.train.psi")||die "can not open";
for($fnum=0;$fnum<@fileList;$fnum++){        
    #print "$fileList[$fnum]\n";
    
    $train_fileList[$fnum] = $fileList[$fnum];    
    @label=();
    @ss1=();    
    @solve=();
    @aa=();
    
    open(IN,"$fileList[$fnum]/$fileList[$fnum].psi")||die "can not open";
    	while($line=<IN>){
    		chomp($line);
    		$line=~s/\s+//g;
				$line=$line/360;
    		push(@label,$line);
    	}
    close IN || die "can not close";
    
    $len1 = @label;
    $trains_lengths[$fnum] = $len1;
    $count = $count + $len1;
    
    open(IN,"$fileList[$fnum]/output1.ss")||die "can not open";  ### ss1
    	while($line=<IN>){
    		chomp($line);
    		$line=~s/^\s+//g;
    		@wds = split(/\s+/,$line);
    		$len=@wds;
    		if($len==6){
    			push(@ss1,"$wds[3]   $wds[4]   $wds[5]");	
    			push(@aa,$wds[1]);
    		}
    	}
    close IN || die "can not close";   
    
    open(IN,"$fileList[$fnum]/seq.exp.sa")||die "can not open";  ### ss1
    	while($line=<IN>){
    		chomp($line);
    		$line=~s/^\s+//g;
    		@wds = split(/\s+/,$line);
    		$len=@wds;
    		if($len==2){
    			push(@solve,"$wds[1]");	
    		}
    	}
    close IN || die "can not close";   
    
    $len = @label;
    $len1=@ss1;
    $len2=@solve;
    if(($len!=$len1)||($len1!=$len2)){
				print "$len1 $len2 $len3\n";
				exit();
    }
    for($i=0;$i<($len);$i++){     
    		  $m = 0;
	        printf OUT "%f ",$label[$i];        
	        for($k=-$window;$k<=$window;$k++){
	            $position=$k+$i;
	            if($position<0||$position>=$len){     #blank residue           	               
	               for($n=1;$n<=4;$n++){	               		
	               		$m++;
	               		printf OUT "%d:0  ",$m;
	               }	              
	               ## spx bank residue
	               for($n=1;$n<=7;$n++){	               		
	               		$m++;
	               		printf OUT "%d:0  ",$m;
	               }	 
	                
	    	    }else{	    	       
		    	    @wds = split(/\s+/,$ss1[$position]);
		    	    for($n=1;$n<=3;$n++){
		    	    	$m++;
		    	    	printf OUT "%d:%f  ",$m,$wds[$n-1];
		    	    }	    	     	   
	    	     	## solvent accessibility
	    	     	$m++;
	    	     	printf OUT "%d:%f  ",$m,$solve[$position];	
	    	     	
	    	     	## spx features
	    	     	my $spx = $spx_hash{$aa[$position]};
              #print "Checking $aa[$position] => spx=$spx\n";
              @wds = split(/\s+/,$spx);
              for($n=1;$n<=7;$n++){
                  $m++;
                  printf OUT "%d:%f  ",$m,$wds[$n-1];
              }
              	    	     	    	     	    	  
	    	    }	    	     
	        }   	
	        print OUT "\n"; 
    }    	    
}
close OUT || die "can not close";
print "real T $count\n";    
    
$count = 0;    
###############################   testing dataset #############################333
@fileList=();
open(IN,"testing_dataset.txt")||die "can not open";
    while($line=<IN>){
        chomp($line);
        push(@fileList,$line);
    }
close IN || die "can not close";
open(OUT,">test_svm.test.psi")||die "can not open";
for($fnum=0;$fnum<@fileList;$fnum++){
    #print "$fileList[$fnum]\n";
	  $test_fileList[$fnum] = $fileList[$fnum];
    @label=();
    @ss1=();
    @solve=();
    @aa=();
    
    open(IN,"$fileList[$fnum]/$fileList[$fnum].psi")||die "can not open";
        while($line=<IN>){
                chomp($line);
                $line=~s/\s+//g;
                $line=$line/360;
                push(@label,$line);
        }
    close IN || die "can not close";
    
    $len1 = @label;
    $tests_lengths[$fnum] = $len1;
    $count = $count + $len1;
    
    open(IN,"$fileList[$fnum]/output1.ss")||die "can not open";  ### ss1
        while($line=<IN>){
                chomp($line);
                $line=~s/^\s+//g;
                @wds = split(/\s+/,$line);
                $len=@wds;
                if($len==6){
                        push(@ss1,"$wds[3]   $wds[4]   $wds[5]");
                        push(@aa,$wds[1]);
                }
        }
    close IN || die "can not close";
    open(IN,"$fileList[$fnum]/seq.exp.sa")||die "can not open";  ### ss1
        while($line=<IN>){
                chomp($line);
                $line=~s/^\s+//g;
                @wds = split(/\s+/,$line);
                $len=@wds;
                if($len==2){
                        push(@solve,"$wds[1]");
                }
        }
    close IN || die "can not close";
    $len = @label;
    $len1=@ss1;
    $len2=@solve;
    if(($len!=$len1)||($len1!=$len2)){
        print "$len1 $len2 $len3\n";
        exit();
    }
    for($i=0;$i<($len);$i++){
                $m = 0;
                printf OUT "%f ",$label[$i];
                for($k=-$window;$k<=$window;$k++){
                    $position=$k+$i;
                    if($position<0||$position>=$len){     #blank residue
                       for($n=1;$n<=4;$n++){
                                $m++;
                                printf OUT "%d:0  ",$m;
                       }
                       
                       ## spx bank residue
					             for($n=1;$n<=7;$n++){	               		
					               		$m++;
					               		printf OUT "%d:0  ",$m;
					             }
					             	 
                    }else{
                            @wds = split(/\s+/,$ss1[$position]);
                            for($n=1;$n<=3;$n++){
                                $m++;
                                printf OUT "%d:%f  ",$m,$wds[$n-1];
                            }
                            ## solvent accessibility
                            $m++;
                            printf OUT "%d:%f  ",$m,$solve[$position];                            
                            
                            ## spx features
								    	     	my $spx = $spx_hash{$aa[$position]};
							              #print "Checking $aa[$position] => spx=$spx\n";
							              @wds = split(/\s+/,$spx);
							              for($n=1;$n<=7;$n++){
							                  $m++;
							                  printf OUT "%d:%f  ",$m,$wds[$n-1];
							              }
							                                          
                    }
                }
                print OUT "\n";
    }
}
close OUT || die "can not close";
print "real TTTT $count\n";
#goto Layer2;

goto L_2;

##############################    end of testing dataset #########################    	
system("./svm-train -s 3 -t 2 train_svm.train.psi my_model");
system("./svm-predict train_svm.train.psi my_model my_output.train");
system("./svm-predict test_svm.test.psi my_model my_output.test");
system("perl MSE_calc.pl train_svm.train.psi my_output.train > MSE.psi.train");
system("perl MSE_calc.pl test_svm.test.psi my_output.test > MSE.psi.test");

test_pos:;

L_2:;

##############################  two layer support verctor machine ##################
my (@data2,@data4);
open(IN,"my_output.train")||die "can not close";
	while($line=<IN>){
			chomp($line);
			push(@data2,$line);
	}
close IN || die "can not close";

open(IN,"my_output.test")||die "can not close";
	while($line=<IN>){
			chomp($line);
			push(@data4,$line);
	}
close IN || die "can not close";

$len2 = @data2;
$len4 = @data4;
print "len2=$len2 len4=$len4\n";
my $temp;
$count = 0;
for($i=0;$i<@train_fileList;$i++){
		$count = $count + $trains_lengths[$i];
}
print "len2 ccc $len2 Check_len=$count\n";
$count = 0;
for($i=0;$i<@test_fileList;$i++){
		$count = $count + $tests_lengths[$i];
}
print "len4 ccc $len4 Check_len=$count\n";
############################# output the predicted value of the first layer  ######
$count = 0;
for($i=0;$i<@train_fileList;$i++){
	  $len1 = $trains_lengths[$i];
		open(OUT,">$train_fileList[$i]/$train_fileList[$i].svm_layer1")||die "can not open";
				for($j=0;$j<$len1;$j++){
						print OUT "$data2[$count]\n";
						$count++;
				}
		close OUT || die "can not close";
}
$count = 0;
for($i=0;$i<@test_fileList;$i++){
	  $len1 = $tests_lengths[$i];
		open(OUT,">$test_fileList[$i]/$test_fileList[$i].svm_layer1")||die "can not open";
				for($j=0;$j<$len1;$j++){
						print OUT "$data4[$count]\n";
						$count++;
				}
		close OUT || die "can not close";
}
############################# end of output values of the first layer ############# 

############################# prepare the secondary layer training adn testing data ##############
@fileList=();
open(IN,"training_dataset.txt")||die "can not open";
    while($line=<IN>){
        chomp($line);
        push(@fileList,$line);
    }
close IN || die "can not close";

open(OUT,">train_svm.train.psi.layer_2")||die "can not open";
for($fnum=0;$fnum<@fileList;$fnum++){        
    #print "$fileList[$fnum]\n";
    
    $train_fileList[$fnum] = $fileList[$fnum];
    
    @label=();
    @ss1=();    
    @solve=();
    @layer1_predict_train=();
    
    open(IN,"$fileList[$fnum]/$fileList[$fnum].psi")||die "can not open";
    	while($line=<IN>){
    		chomp($line);
    		$line=~s/\s+//g;
				$line=$line/360;
    		push(@label,$line);
    	}
    close IN || die "can not close";
    
    $len1 = @label;
    $trains_lengths[$fnum] = $len1;
    
    open(IN,"$fileList[$fnum]/output1.ss")||die "can not open";  ### ss1
    	while($line=<IN>){
    		chomp($line);
    		$line=~s/^\s+//g;
    		@wds = split(/\s+/,$line);
    		$len=@wds;
    		if($len==6){
    			push(@ss1,"$wds[3]   $wds[4]   $wds[5]");	
    		}
    	}
    close IN || die "can not close";   
    
    open(IN,"$fileList[$fnum]/seq.exp.sa")||die "can not open";  ### ss1
    	while($line=<IN>){
    		chomp($line);
    		$line=~s/^\s+//g;
    		@wds = split(/\s+/,$line);
    		$len=@wds;
    		if($len==2){
    			push(@solve,"$wds[1]");	
    		}
    	}
    close IN || die "can not close";   
    
    open(IN,"$fileList[$fnum]/$fileList[$fnum].svm_layer1")||die "can not open";  ### ss1
    	while($line=<IN>){
    		push(@layer1_predict_train,$line);
    	}
    close IN || die "can not close";   
    
    $len = @label;
    $len1=@ss1;
    $len2=@solve;
    $len3=@layer1_predict_train;
    
    if(($len!=$len1)||($len1!=$len2)||($len2!=$len3)){
				print "$len $len1 $len2 $len3\n";
				exit();
    }
    for($i=0;$i<($len);$i++){     
    			$m = 0;
	        printf OUT "%f ",$label[$i];        
	        for($k=-$window;$k<=$window;$k++){
	            $position=$k+$i;
	            if($position<0||$position>=$len){     #blank residue           	               
	               for($n=1;$n<=5;$n++){	               		
	               		$m++;
	               		printf OUT "%d:0  ",$m;
	               }	               
	    	    }else{	    	       
		    	    @wds = split(/\s+/,$ss1[$position]);
		    	    for($n=1;$n<=3;$n++){
		    	    	$m++;
		    	    	printf OUT "%d:%f  ",$m,$wds[$n-1];
		    	    }	    	     	   
	    	     	## solvent accessibility
	    	     	$m++;
	    	     	printf OUT "%d:%f  ",$m,$solve[$position];	 
	    	     	
	    	     	## predicted value
	    	     	$m++;
	    	     	printf OUT "%d:%f  ",$m,$layer1_predict_train[$position];   	     	    	  
	    	    }	    	     
	        }   	
	        print OUT "\n"; 
    }    	    
}
close OUT || die "can not close";
    
###############################   testing dataset #############################333
@fileList=();
open(IN,"testing_dataset.txt")||die "can not open";
    while($line=<IN>){
        chomp($line);
        push(@fileList,$line);
    }
close IN || die "can not close";

open(OUT,">test_svm.test.psi.layer_2")||die "can not open";
for($fnum=0;$fnum<@fileList;$fnum++){
    #print "$fileList[$fnum]\n";
	  $test_fileList[$fnum] = $fileList[$fnum];
    @label=();
    @ss1=();
    @solve=();
    @layer1_predict_test=();
    
    open(IN,"$fileList[$fnum]/$fileList[$fnum].psi")||die "can not open";
        while($line=<IN>){
                chomp($line);
                $line=~s/\s+//g;
                $line=$line/360;
                push(@label,$line);
        }
    close IN || die "can not close";
    
    $len1 = @label;
    $tests_lengths[$fnum] = $len1;
    
    open(IN,"$fileList[$fnum]/output1.ss")||die "can not open";  ### ss1
        while($line=<IN>){
                chomp($line);
                $line=~s/^\s+//g;
                @wds = split(/\s+/,$line);
                $len=@wds;
                if($len==6){
                        push(@ss1,"$wds[3]   $wds[4]   $wds[5]");
                }
        }
    close IN || die "can not close";
    
    open(IN,"$fileList[$fnum]/seq.exp.sa")||die "can not open";  ### ss1
        while($line=<IN>){
                chomp($line);
                $line=~s/^\s+//g;
                @wds = split(/\s+/,$line);
                $len=@wds;
                if($len==2){
                        push(@solve,"$wds[1]");
                }
        }
    close IN || die "can not close";
    
    open(IN,"$fileList[$fnum]/$fileList[$fnum].svm_layer1")||die "can not open";  ### ss1
    	while($line=<IN>){
    		push(@layer1_predict_test,$line);
    	}
    close IN || die "can not close"; 
    
    $len = @label;
    $len1=@ss1;
    $len2=@solve;
    $len3=@layer1_predict_test;
    
    if(($len!=$len1)||($len1!=$len2)||($len2!=$len3)){
        print "$len $len1 $len2 $len3\n";
        exit();
    }
    for($i=0;$i<($len);$i++){
                $m = 0;
                printf OUT "%f ",$label[$i];
                for($k=-$window;$k<=$window;$k++){
                    $position=$k+$i;
                    if($position<0||$position>=$len){     #blank residue
                       for($n=1;$n<=5;$n++){
                                $m++;
                                printf OUT "%d:0  ",$m;
                       }
                    }else{
                            @wds = split(/\s+/,$ss1[$position]);
                            for($n=1;$n<=3;$n++){
                                $m++;
                                printf OUT "%d:%f  ",$m,$wds[$n-1];
                            }
                            ## solvent accessibility
                            $m++;
                            printf OUT "%d:%f  ",$m,$solve[$position];
                            $m++;
                            printf OUT "%d:%f  ",$m,$layer1_predict_test[$position];
                    }
                }
                print OUT "\n";
    }
}
close OUT || die "can not close";
############################# end of preparing data ##############################################
system("./svm-train -s 3 -t 2 train_svm.train.psi.layer_2 my_model_layer2");
system("./svm-predict test_svm.test.psi.layer_2 my_model_layer2 my_output.test.layer2");
system("perl MSE_calc.pl test_svm.test.psi my_output.test.layer2 > lay2.result");

