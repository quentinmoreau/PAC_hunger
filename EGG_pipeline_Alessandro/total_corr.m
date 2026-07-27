function TC=total_corr(x1, x2, x3, x4, a)

%trasformo tutti in vettori riga se non lo sono
if size(x1, 1)>size(x1, 2)
    x1=x1';
end
if size(x2, 1)>size(x2, 2)
    x2=x2';
end
if size(x3, 1)>size(x3, 2)
    x3=x3';
end
if size(x4, 1)>size(x4, 2)
    x4=x4';
end

%Silverman's Rule of Thumb
s1=1.06*std(x1);
%{
%Gram matrix for x1
G1=zeros(length(x1), length(x1));
for i=1:size(G1, 1)
    for j=1:size(G1, 2)
        G1(i, j)= exp(-( x1(i)-x1(j) )^2/(2*s1) );
    end
end
%}
dist_sq = bsxfun(@minus, x1, x1').^2;
G1 = exp(-dist_sq / (2*s1));

%Silverman's Rule of Thumb
s2=1.06*std(x2);
%{
%Gram matrix for x1
G2=zeros(length(x2), length(x2));
for i=1:size(G2, 1)
    for j=1:size(G2, 2)
        G2(i, j)= exp(-( x2(i)-x2(j) )^2/(2*s2) );
    end
end
%}
dist_sq = bsxfun(@minus, x2, x2').^2;
G2 = exp(-dist_sq / (2*s2));

%Silverman's Rule of Thumb
s3=1.06*std(x3);
%{
%Gram matrix for x1
G3=zeros(length(x3), length(x3));
for i=1:size(G3, 1)
    for j=1:size(G3, 2)
        G3(i, j)= exp(-( x3(i)-x3(j) )^2/(2*s3) );
    end
end
%}

dist_sq = bsxfun(@minus, x3, x3').^2;
G3 = exp(-dist_sq / (2*s3));



%Silverman's Rule of Thumb
s4=1.06*std(x4);

%{
%Gram matrix for x1
G4=zeros(length(x4), length(x4));
for i=1:size(G4, 1)
    for j=1:size(G4, 2)
        G4(i, j)= exp(-( x4(i)-x4(j) )^2/(2*s4) );
    end
end
%}

dist_sq = bsxfun(@minus, x4, x4').^2;
G4 = exp(-dist_sq / (2*s4));

%{
%A matrix for series x1
A1=zeros(size(G1));
for i=1:size(A1, 1)
    for j=1:size(A1, 2)
        A1(i, j)= G1(i, j)/sqrt(G1(i, i)*G1(j, j));
    end
end
A1=A1/(size(G1, 1));
%}

A1 = zeros(size(G1));
diagG1 = diag(G1);

% Operazioni vettorizzate
A1 = bsxfun(@rdivide, G1, sqrt(diagG1'*diagG1));

% Normalizzazione per il numero di righe
A1 = A1 / size(G1, 1);

%{
%A matrix for series x2
A2=zeros(size(G2));
for i=1:size(A2, 1)
    for j=1:size(A2, 2)
        A2(i, j)= G2(i, j)/sqrt(G2(i, i)*G2(j, j));
    end
end
A2=A2/(size(G2, 1));
%}

A2 = zeros(size(G2));
diagG2 = diag(G2);

% Operazioni vettorizzate
A2 = bsxfun(@rdivide, G2, sqrt(diagG2'*diagG2));

% Normalizzazione per il numero di righe
A2 = A2 / size(G2, 1);


%{
%A matrix for series x3
A3=zeros(size(G3));
for i=1:size(A3, 1)
    for j=1:size(A3, 2)
        A3(i, j)= G3(i, j)/sqrt(G3(i, i)*G3(j, j));
    end
end
A3=A3/(size(G3, 1));
%}


A3 = zeros(size(G3));
diagG3 = diag(G3);

% Operazioni vettorizzate
A3 = bsxfun(@rdivide, G3, sqrt(diagG3'*diagG3));

% Normalizzazione per il numero di righe
A3 = A3 / size(G3, 1);


%{
%A matrix for series x4
A4=zeros(size(G4));
for i=1:size(A4, 1)
    for j=1:size(A4, 2)
        A4(i, j)= G4(i, j)/sqrt(G4(i, i)*G4(j, j));
    end
end
A4=A4/(size(G4, 1));
%}

A4 = zeros(size(G4));
diagG4 = diag(G4);

% Operazioni vettorizzate
A4 = bsxfun(@rdivide, G4, sqrt(diagG4'*diagG4));

% Normalizzazione per il numero di righe
A4 = A4 / size(G4, 1);
%G4

%trace(A4)
%trace(A4^2)

if a==2
    tA1=A1(:).'*reshape(A1.',[], 1);
    tA2=A2(:).'*reshape(A2.',[], 1);
    tA3=A3(:).'*reshape(A3.',[], 1);
    tA4=A4(:).'*reshape(A4.',[], 1);
    p12=A1*A2;
    p34=A3*A4;
    t1234=p12(:).'*reshape(p34.',[], 1);
    p1234=A1*A2*A3*A4;
    t1234_2=p1234(:).'*reshape(p1234.',[], 1);

    TC=(1/(1-a))*log2(tA1)+(1/(1-a))*log2(tA2)+...
    (1/(1-a))*log2(tA3)+(1/(1-a))*log2(tA4)-...
    (1/(1-a))*log2(t1234_2/(t1234^2));

else

TC=(1/(1-a))*log2(trace(A1^a))+(1/(1-a))*log2(trace(A2^a))+...
    (1/(1-a))*log2(trace(A3^a))+(1/(1-a))*log2(trace(A4^a))-...
    (1/(1-a))*log2(trace( (A1*A2*A3*A4/(trace(A1*A2*A3*A4)))^a ));
end



end