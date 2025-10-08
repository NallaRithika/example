#include<stdio.h>
#include<conio.h>
#define SIZE 10
void push(int);
void pop();
void display();
int stack[SIZE],top=-1;
void main()
{
   int value,choice;
   clrscr();
   while(1)
   {
       printf("\n\n**menu**\n\n");
       printf("1.push\n 2.pop\n 3.display\n 4.exit");
       printf("\n enter your choice");
       scanf("%d",choice);
       switch(choice)
       {
           case 1:printf("enter value to be inserted");
           scanf("%d",&value);
           push(value);
           break;
           case 2:pop();
           break;
           case 3:display();
           break;
           case 4:exit(0);
           default:printf("\n wrong seletion");
           }
   }
}
void push(int value)
{
    if(top==SIZE-1)
       printf("\n stack is full");
   else
    {
       top++;
       stack[top]=value;
       printf("\n insertion success");
   }
}
void pop()
{
if(top==-1)
   printf("\n stack empty deletion not possible");
else
{
   printf("\n deleted %d",stack[top]);
   top--;
   }
}
void display()
{
    if(top==-1)
        printf("\n stack empty");
    else
    {
        int i;
        printf("stack element are");
        for(i=top;i>=0;i--)
            printf("%d\n",stack[i]);
        }
}
