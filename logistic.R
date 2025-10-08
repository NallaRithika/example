#sample data
x<-c(1,2,3,4,5,1,2,3,4,5)
y<-c(0,0,0,1,1,0,1,1,1,1)
#create a data frame
df<-data.frame(x,y)
model<-glm(y~x,family=binomial(link="logit"),data=df)
summary(model)
new_x<-data.frame(x=c(2.5,3.5))
predicted_probabilities<-predict(model,newdata=new_x,type="response")
print(predicted_probabilities)
predicted_classes<-ifelse(predicted_probabilities>0.5,"1","0")
print(predicted_classes)

