

x<-c(1,2,3,4,5)
y<-c(2,4,5,4,5)
df<-data.frame(x,y)
model<-lm(y~x,data=df)
summary(model)
new_x<-data.frame(x=6)
predicted_y<-predict(model,newdata=new_x)
print(predicted_y)
plot(x,y,main="LinearRegression",xlab="x-axis",ylab="y-axis")
abline(model,col="red")
