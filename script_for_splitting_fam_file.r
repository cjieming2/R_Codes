let i=0 a=14
     let j=0

     while [ $i -le $a ]
     do
      while [ $j -le $a ]
      do

         bsub -o /dev/null -e /dev/null plink --bfile plink\
                 --all \
                 --genome \
                 --genome-lists tmp.list`printf "%03i\n" $i` \
                                tmp.list`printf "%03i\n" $j` \
                 --out data.sub.$i.$j

      let j=$j+1
      done

      let i=$i+1
      let j=$i

     done
