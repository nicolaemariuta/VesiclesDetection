function Icorrected = correctSlice( I,n )
%I - input 3D image
%n - slice to repair
Icorrected = I;

%make correction
a=sum(sum(Icorrected(:,:,n-1)));
b=sum(sum(Icorrected(:,:,n+1)));
c=sum(sum(Icorrected(:,:,n)));
Icorrected(:,:,37)=Icorrected(:,:,37)*(a+b)/2/c;






