function [FFDistance] = FFDistance(Molecule1Coords,row1,row2,row3,row4)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pos1=Molecule1Coords(row1,:);
pos2=Molecule1Coords(row2,:);
pos3=Molecule1Coords(row3,:);
pos4=Molecule1Coords(row4,:);

FFDistance13=(((pos1(1,1)-pos3(1,1))^2)+((pos1(1,2)-pos3(1,2))^2)+((pos1(1,3)-pos3(1,3))^2))^(1/2)
FFDistance14=(((pos1(1,1)-pos4(1,1))^2)+((pos1(1,2)-pos4(1,2))^2)+((pos1(1,3)-pos4(1,3))^2))^(1/2)
FFDistance23=(((pos2(1,1)-pos3(1,1))^2)+((pos2(1,2)-pos3(1,2))^2)+((pos2(1,3)-pos3(1,3))^2))^(1/2)
FFDistance24=(((pos2(1,1)-pos4(1,1))^2)+((pos2(1,2)-pos4(1,2))^2)+((pos2(1,3)-pos4(1,3))^2))^(1/2)

end

