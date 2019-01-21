# Input physlog file. Make sure files follow naming convention: fileprefix#.log where # is the participant ID number
# Output cardiac position file, flag file, and resp file. Heartrate and total volumes output to command line (HR in 1st col, TV in 2nd col)
#makephysiofiles <- function(num,volumes=230,slices=28,Directory,fileprefix)
#Setwd as the foldetotr with all your physlog files.
#Athena Theyers modified by Sarah Atwi Tuesday, January 8th, 2019

num=c()
volumes =  # Number of volumes in fMRI scan 
slices =  # Number of slices in fMRI scan
heartrate = vector(length=length(num))
totalvolumes = vector(length=length(num))#check if it meets your thresholded. 
allflags = matrix(0,nrow=volumes,ncol=length(num))
InputDirectory = paste0("", num, '.log') #source directory of where physlogs are
fileprefix = '' #file prefix name if you have it. leave blank.
OutputDirectory = paste0("")

for (m in 1:length(num))
{
	name = paste0("", num[m],'.log') # physlog file
	filedata = read.table(name,skip=5,sep="")
	ppu = filedata[,5]
	resp = filedata[,6]
	gx = filedata[,7]
	mark = filedata[,10]
	cardiacposition <- vector(length=length(mark))

	point <- 1
	maxlength <-vector(length=0)
	for (i in 1:length(cardiacposition))
	{
		if (mark[i]==2|mark[i]==22) #find marker for end of each cardiac cycle
		{
			maxlength <- c(maxlength,point)
			point <- 0
		}else{
			point <- point +1
		}
		cardiacposition[i] <- point
	}

	heads <- 1
	cardiacposition1 <- vector(length=length(cardiacposition))

	for (i in 1:(length(cardiacposition)-1))
	{
		tails <- i
		if (cardiacposition[tails+1]==0)
		{
			for(j in heads:tails)
			{
				cardiacposition1[j] <- cardiacposition[j]/cardiacposition[tails]
			}
			heads <- tails + 2
		}
	}

	tails <- tails + 1

	for(j in heads:tails)
		{
			cardiacposition1[j] <- cardiacposition[j]/cardiacposition[tails]
		}

	ind2 <- vector(length=0)
	i <- 2
	while(i<length(gx)-1)
	{
		if (gx[i]> 1500& gx[i-1]<=gx[i]&gx[i+1]<=gx[i]) #find readout of slice. If missing slices, reduce gx limit, default: gx[i] > 1500
		{
			ind2 <- c(ind2,i)
			i <- i+1		
		}
		i <- i+1
	}
	ind3 <- vector(length=slices*volumes)

	for (i in 1:volumes)
	{
		for (j in 1:slices)
		{
			ind3[slices*(i-1)+j] <- ind2[length(ind2)-(slices-j)-slices*(volumes-i)] 
		}
	}
	namecard = paste0(OutputDirectory,fileprefix,"Cardiac",num[m],".txt") # Cardiac position file
#write(cardiacposition1[ind3],file=namecard,ncolumns=slices)
	nameresp = paste0(OutputDirectory,fileprefix,"Resp",num[m],".txt") # Resp file
#write(resp[ind3[1]:length(resp)],file=nameresp,ncolumns=1)
	flags = matrix(0,nrow=volumes,ncol=slices)
	temp = cardiacposition[ind3]/cardiacposition1[ind3]
	temp[is.na(temp)] = 0
	standard = mean(temp[temp<700&temp>100])
	deviation = sd(temp[temp<700&temp>100])
	temp = matrix((temp),nrow=slices,ncol=volumes)
  temp=t(temp)
	if (2*deviation > 150)
	{
		upperlimit = standard + 2*deviation
		lowerlimit = standard - 2*deviation
	}else{
		upperlimit = standard + 150
		lowerlimit = standard - 150
	}
 
	flags = temp<upperlimit&temp>lowerlimit|temp==0
	flags[is.na(flags)] = 1
	allflags[,m] = flags[,1]&flags[,2]&flags[,3]&flags[,4]&flags[,5]&flags[,6]&flags[,7]&flags[,8]&flags[,9]&
    flags[,10]&flags[,11]&flags[,12]&flags[,13]&flags[,14]&flags[,15]&flags[,16]&flags[,17]&flags[,18]&
    flags[,19]&flags[,20]&flags[,21]&flags[,22]&flags[,23]&flags[,24]&flags[,25]&flags[,26]&flags[,27]&flags[,28]
  totalvolumes[m] = sum(allflags[,m]==1)
	heartrate[m] = standard/500*60 # Calculate average heartrate in bpm
message(num[m])
}

#view total volumes
View(t(totalvolumes))

#write the flagged volumes
nameflags = paste0(".txt") # Flag file
write(t(allflags),file=nameflags,ncolumns=length(num))

###### to remove flags
Remove = c() # check out the output flag file and see the participants to remove because of poor cardiac. 
#Numbers should be order in which participants were run NOT the participant ID numbers. e.g. ID = 5,7,24,32,57, to remove 7 and 32, Remove = c(2,4)
flags <- read.table(nameflags,header=FALSE,sep="")
flags_removed <- subset(flags,select=-Remove)
write.table(flags_removed,paste0(".txt"),row.name=FALSE,col.name=FALSE)
