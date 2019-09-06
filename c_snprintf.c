#include<stdio.h>
#include<string.h>
int main()
{
    char s[]="XXXXXXXXXXXXXXXXXXXX";
    printf("Original string=%s\t",s);
    printf("Length=%d\n",strlen(s));
    int cx=snprintf(s,10,"HelloWorldIsTruncted");
    printf("cx=%d\n",cx);
    printf("Truncated string=%s\tLength=%d\n",s,strlen(s));
    s[9]='Y';
    s[20]='\0';
    printf("Ruined string=%s\n",s);
}

