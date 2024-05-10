cd $1 && for i in `ls | sort -V`; 
do 
	ln -s $1/$i $2/$i; 
	sleep 5 # sleep for 5 seconds
done
