betalabel<-c("crent",
	"distance 10 min vs 5 min",
	"distance 15 min vs 5 min",
	"noise	moderate vs bad",
	"noise	very good vs bad",
	"facilities private bath & toilet, shared kitchen vs all shared",
	"facilities private toilet, shared shower & kitchen vs all shared",
	"facilities private toilet, shower & kitchen vs all shared",
	"condition newly renovated vs not renovated")


        
label<-c("rent 250",
	"rent 350",
	"rent 450",
	"distance 5 min",
	"distance 10 min",
	"distance 15 min",
	"noise 	bad",
	"noise	moderate",
	"noise	very good",
	"facilities shared toilet, shower & kitchen",
	"facilities private bath & toilet, shared kitchen",
	"facilities private toilet, shared shower & kitchen",
	"facilities private toilet, shower & kitchen",
	"condition not renovated",
	"condition newly renovated")

#consider
I<-107
M<-5
S<-14

consider<-matrix(scan("C:/Users/u0105757/Desktop/code/subset conjunctive model/data/consider_student.txt"),nrow=107,byrow=T)

cons<-array(consider,c(I,M,S))
#apply(cons,c(1,3),sum)

#read choice data
cstudent<-matrix(scan("C:/Users/u0105757/Desktop/code/subset conjunctive model/data/choice_student.txt"),nrow=107,byrow=T)
# drop student ID
cstudent<-cstudent[,-1]

# create array of binary choice data
#last category is opt-out
I<-107
S<-14
N<-6
student<-array(rep(0,I*S*N),c(I,S,N))

for (i in 1:I){
 for (s in 1:S){
  if (cstudent[i,s]==1) {student[i,s,1]<-1}
  else if (cstudent[i,s]==2) {student[i,s,2]<-1}
  else if (cstudent[i,s]==3) {student[i,s,3]<-1}
  else if (cstudent[i,s]==4) {student[i,s,4]<-1}
  else if (cstudent[i,s]==5) {student[i,s,5]<-1}
}}


# read design
# matrix of 70 choice sets by 15 attribute levels
# order of choice sets: alternative 1 set 1, alternative 2 set 1, ..., alt 5 set 1, alt 1 set 2, ...
# order of attribute levels: see vector of attribute labels "label" 

des<-matrix(scan("C:/Users/u0105757/Desktop/code/subset conjunctive model/data/design_binair.txt"),nrow=70,byrow=T)

#price in 100 euro
price<-des[,1:3]%*%c(2.5, 3.5, 4.5)

# select design for regression model
# include continuous price
# include per attribute K_p-1 attribute levels
seldesign<-cbind(price,des[,c(5,6,8,9,11,12,13,15)])


# write.table(seldesign,file="C:\\data_HUB\\projecten\\89_subset_conjunctive_choice\\datasets\\seldesign.txt",sep=" ",row.names=FALSE,col.names=FALSE)

