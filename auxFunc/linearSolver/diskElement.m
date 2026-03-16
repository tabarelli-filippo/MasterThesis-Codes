function [M0e,C1e] = diskElement(Disk)
%DISKELEMENT computes the element matrix for a given disk description.
% For each element see Manual
%
%INPUT: Disk   Structure
%
%OUTPUT:M0e  mass matrix
%       C1e  gyroscopic matrix

type = Disk.type;
switch type
    case 1
        D_ext = Disk.D_ext;
        D_int = Disk.D_int;
        thick = Disk.thick;
        rho = Disk.rho;
        
        m = rho * pi / 4 * thick * (D_ext^2 - D_int^2);
        Ip = rho * pi / 32 * thick * (D_ext^4 - D_int^4);
        Id = Ip / 2 + m * (thick^2) / 12;
    case 2
        m = Disk.mass;
        Ip = Disk.Ip;
        Id = Disk.Id;
end
    M0e = diag([m m Id Id]);
    C1e = zeros(4,4);
    C1e(3,4) = Ip;
    C1e(4,3) = -Ip;
end