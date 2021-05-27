function [Am, Bm, Cm, Dm, Al, Bl, Cl, Dl, length] = e223115_aktas(text_file,input_file)
%% PART-1

lib = readtable(input_file);                   
conductor_names = lib.CodeWord;               

fid = fopen(text_file);                          

if fid == -1                        
    error("The file could not be found.");      
    return
end

tline = fgetl(fid);                
raw_data = convertCharsToStrings(tline);
while ischar(tline)                               
    tline = fgetl(fid);                          
    if strcmp(tline , "-999")
        break
    end
    raw_data=[raw_data ; convertCharsToStrings(tline)];
end
fclose(fid);
    while ~any(strcmp(conductor_names,raw_data(10)))    
        warning('The conductor does not exist. Try another conductor name.');
        raw_data(10) = input('Conductor name: ','s');            
    end                                          
n_cct = str2double(raw_data(2));                          
n_bundle = str2double(raw_data(4));                        
d_bundle = str2double(raw_data(6));
length = str2double(raw_data(8));
line_type = raw_data(10);
   
%% PART-2
iRow = (strcmp(conductor_names, line_type)==1);
lib_datas = lib(iRow,:);
r = table2array(lib_datas(1,8)) * 0.3048; %Line specific CMR in terms of m
rc= table2array(lib_datas(1,5)) * 0.0254 * 0.5;%convert to m
r_b = r;
r_c = rc;

coords1 = [(str2double(raw_data(12)) + 1i*str2double(raw_data(13)))...
           (str2double(raw_data(15)) + 1i*str2double(raw_data(16)))...
           (str2double(raw_data(18)) + 1i*str2double(raw_data(19)))];
coords1_images = conj(coords1);

if d_bundle ~= 0
    angle_bundle = pi/n_bundle;
    R = d_bundle/(2*sin(angle_bundle));
    % shifted version of bundle coords. Centered at 0,0
    bundle_coords = R*exp( 1i * (0:n_bundle-1)' * 2 * angle_bundle);
    r_b = r_b * prod(  abs(bundle_coords(1) - bundle_coords(2:end))  );
    r_c = rc * prod(  abs(bundle_coords(1) - bundle_coords(2:end))  );
    r_b = (r_b)^(1/n_bundle); 
    r_c = (r_c)^(1/n_bundle);    
else
    bundle_coords = 0;
end
        all_bundle_coords = [(bundle_coords + coords1(1))...
    (bundle_coords + coords1(2))     bundle_coords + coords1(3)];

        bundle_coords1_images = [(bundle_coords + coords1_images(1)) ...
(bundle_coords + coords1_images(2))  ( bundle_coords + coords1_images(3))];

    A=[repmat(conj(all_bundle_coords(:,1))',n_bundle,1);...
    repmat(conj(all_bundle_coords(:,2))',n_bundle,1);...
    repmat(conj(all_bundle_coords(:,3))',n_bundle,1)];
    [m, n] = size(A);

d_ab = abs(coords1(1)-coords1(2));
d_ac = abs(coords1(1)-coords1(3));
d_bc = abs(coords1(2)-coords1(3));

if n_cct == 1
    % GMD Bonus calculation
    B=[all_bundle_coords(:,2);...
    all_bundle_coords(:,3);...
    all_bundle_coords(:,1)];

    GMD_b=prod(prod(abs(A-B)).^(1/max(m,n)))^(1/min(m,n));

    % H_ij Bonus calculation
    B=[bundle_coords1_images(:,2);...
    bundle_coords1_images(:,3);...
    bundle_coords1_images(:,1)];

    H_ijb=prod(prod(abs(A-B)).^(1/max(m,n)))^(1/min(m,n));

    % H_i Bonus calculation
    B=[bundle_coords1_images(:,1);...
    bundle_coords1_images(:,2);...
    bundle_coords1_images(:,3)];

    H_ib=prod(prod(abs(A-B)).^(1/max(m,n)))^(1/min(m,n));

    % image distance calculations:
    H_i  = [abs(coords1(1) - coords1_images(1)) ...
abs(coords1(2) - coords1_images(2))   abs(coords1(3) - coords1_images(3))];

    H_ij = [abs(coords1(1) - coords1_images(2)) ...
abs(coords1(3) - coords1_images(1))   abs(coords1(2) - coords1_images(3))];
% effect of earth on capacitance:

    H_n = ((prod(H_ij)/prod(H_i))^(1/3));
    H_b = ((H_ijb)/(H_ib));
    GMD_n = (  d_ab * d_ac * d_bc  )^(1/3);
    GMR = r_b;
    Dsc=r_c;
else
    clear i j;
    
    coords2 = [(str2double(raw_data(21)) + 1i*str2double(raw_data(22)))...
               (str2double(raw_data(24)) + 1i*str2double(raw_data(25)))...
               (str2double(raw_data(27)) + 1i*str2double(raw_data(28)))];
    coords2_images = conj(coords2);
    bundle_coords2_images = bundle_coords + coords2_images;
    
    all_bundle_coords2 = [(bundle_coords + coords2(1))...
    (bundle_coords + coords2(2))     bundle_coords + coords2(3)];
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %A --> [a;a;a';a';b;b;b';b';c;c;c';c'](Transpozed)
    A=[repmat(conj(all_bundle_coords(:,1))',2*n_bundle,1);...
       repmat(conj(all_bundle_coords2(:,1))',2*n_bundle,1);...
       repmat(conj(all_bundle_coords(:,2))',2*n_bundle,1);...
       repmat(conj(all_bundle_coords2(:,2))',2*n_bundle,1);...
       repmat(conj(all_bundle_coords(:,3))',2*n_bundle,1);...
       repmat(conj(all_bundle_coords2(:,3))',2*n_bundle,1)];
   [m, n] = size(A);
    % GMD_b calculation:
    %B --> [b;b';b;b';  c;c';c;c';  a;a';a;a']
    B=[all_bundle_coords(:,2);...
       all_bundle_coords2(:,2);...
       all_bundle_coords(:,2);...
       all_bundle_coords2(:,2);...
    
       all_bundle_coords(:,3);...
       all_bundle_coords2(:,3);...
       all_bundle_coords(:,3);...
       all_bundle_coords2(:,3);...
       
       all_bundle_coords(:,1);...
       all_bundle_coords2(:,1);...
       all_bundle_coords(:,1);...
       all_bundle_coords2(:,1)];
   
   GMD_b = prod(prod(abs(A-B)).^(1/max(m,n)))^(1/min(m,n)); 
    % H_ij bonus calculation    
    B=[bundle_coords1_images(:,2);...
       bundle_coords2_images(:,2);...
       bundle_coords1_images(:,2);...
       bundle_coords2_images(:,2);...
       
       bundle_coords1_images(:,3);...
       bundle_coords2_images(:,3);...
       bundle_coords1_images(:,3);...
       bundle_coords2_images(:,3);...
       
       bundle_coords1_images(:,1);...
       bundle_coords2_images(:,1);...
       bundle_coords1_images(:,1);...
       bundle_coords2_images(:,1)];
    
    H_ijb = prod(prod(abs(A-B)).^(1/max(m,n)))^(1/min(m,n));
    % H_i bonus calculation  
    B=[bundle_coords1_images(:,1);...
       bundle_coords2_images(:,1);...
       bundle_coords1_images(:,1);...
       bundle_coords2_images(:,1);...
    
       bundle_coords1_images(:,2);...
       bundle_coords2_images(:,2);...
       bundle_coords1_images(:,2);...
       bundle_coords2_images(:,2);...
       
       bundle_coords1_images(:,3);...
       bundle_coords2_images(:,3);...
       bundle_coords1_images(:,3);...
       bundle_coords2_images(:,3)];
    
    H_ib = prod(prod(abs(A-B)).^(1/max(m,n)))^(1/min(m,n));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % H Calculation normal part 
        H_i  = [abs(coords1(1) - coords1_images(1)) ...
                   abs(coords1(2) - coords1_images(2)) ...
                   abs(coords1(3) - coords1_images(3))];
        H_ij = [abs(coords1(1) - coords1_images(2))*...
abs(coords1(1) - coords2_images(2))*abs(coords2(1) - coords1_images(2))*...
                 abs(coords2(1) - coords2_images(2))
                 
abs(coords1(1) - coords1_images(3))*abs(coords1(1) - coords2_images(3))*...
abs(coords2(1) - coords1_images(3))*abs(coords2(1) - coords2_images(3))
                 
abs(coords1(2) - coords1_images(3))*abs(coords1(2) - coords2_images(3))*...
abs(coords2(2) - coords1_images(3))*abs(coords2(2) - coords2_images(3))];
  
D = [abs(coords1(1)- coords2_images(1)) abs(coords1(2)- coords2_images(2))...
        abs(coords1(3)- coords2_images(3))];
        H_i= prod(H_i)^(1/3);
        H_i =(prod(sqrt(D*H_i)))^(1/3);
        H_ij= prod(H_ij)^(1/12);         
        H_n = (H_ij/H_i);
        H_b = (H_ijb/H_ib); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    d_ab = abs(coords1(1)-coords1(2))*abs(coords1(1)-coords2(2))*...
        abs(coords2(1)-coords1(2))*abs(coords2(1)-coords2(2));
    d_ac = abs(coords1(1)-coords1(3))*abs(coords1(1)-coords2(3))*...
        abs(coords2(1)-coords1(3))*abs(coords2(1)-coords2(3));
    d_bc = abs(coords1(2)-coords1(3))*abs(coords1(2)-coords2(3))*...
        abs(coords2(2)-coords1(3))*abs(coords2(2)-coords2(3));
    
    D = [abs(coords1(1)-coords2(1)) abs(coords1(2)-coords2(2))...
        abs(coords1(3)-coords2(3))];
    GMR = (prod(sqrt(D*r_b)))^(1/3);
    Dsc = (prod(sqrt(D*r_c)))^(1/3);
    GMD_n = (d_ab*d_ac*d_bc)^(1/12);
    % dist ab
    all_bundle_coords(:,1);
    repmat(all_bundle_coords(:,1),1,n_bundle);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R        =  table2array(lib_datas(1,7))*0.62137/(n_bundle*n_cct); %ohm/ml to ohm/km
X        =  (2*pi*50)*2*10^(-7+3)*log(GMD_n/GMR);
X_Bonus  =  (2*pi*50)*2*10^(-7+3)*log(GMD_b/GMR);
B        =  (2*pi*50)*(2*pi*8.854*10^(-12+3))/(log(GMD_n/Dsc)-log(H_n));
B_Bonus  =  (2*pi*50)*(2*pi*8.854*10^(-12+3))/(log(GMD_b/Dsc)-log(H_b));

%% PART-3

Xm = 1i*X;
Ym = 1i*B;
Z = (R + Xm) * length;
Y = Ym * length;

Am = 1 + (Z * Y) / 2;
Bm = Z;
Cm = Y * (1 + (Z * Y) / 4);
Dm = Am;

z0 = sqrt(Z/Y);
al = sqrt(Z*Y);
Al = cosh(al);
Bl = z0*sinh(al);
Cl = sinh(al)/z0;
Dl = Al;
end
