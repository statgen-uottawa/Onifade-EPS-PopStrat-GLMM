#!/bin/bash

for i in {1..1000}; do
	
	awk '{ if ( $3 == 2 ) { $3 = 1 } else if ( $3 == 1 ) { $3 = 0 }; print}' EPS_pheno${i}.phe > pheno${i}.phe

done


