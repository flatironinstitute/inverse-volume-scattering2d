function [sensor_info,data] = C2Mtranslator(sensors,directions,u_meas)
% This code translates the structs from Carlos code 2 Manas code for single
% frequency data and information of sensors and directions
% INPUT
%      sensors -> struct of the form sensors.direction(id).position, where id is the
%      direction id in the vector directions and sensors(id).position is a
%      2xNr_{id} vector. Note that the number of receptors can change for each
%      direction
%      directions -> 2xNd vector
%      umeas -> struct of the form umeas.direction(id).field, where id is
%      the direction direction(:,id) in the vector directions and
%      umeas.direction(id).field is a complex vector with the measuemnt of
%      the data at the equivalent sensor from sensors(id).position
%
% OUTPUT
% sensor_info - sensor information struct
%               sensor_info.tgt(2,nmeas) - xy cooordinates of sensors
%               sensor_info.tgt(1:2,i) = xy coordinates corresponding the 
%               ith measurement
%               sensor_info.t_dir(nmeas) - incident directions
%               sensor_info.t_dir(i) - is the incident direction 
%               corresponding to the ith measurement
% data - vector with data for the sensors in sensor_info.tgt pretty much
%        creates a long vector with all the data

ntheta = size(directions,2); %number of directions

sensor_info.tgt = [];
sensor_info.t_dir = [];
data = [];
for id = 1 : ntheta
    
    nsensors = size(sensors.direction(id).position,2);
    sensor_info.tgt = [sensor_info.tgt sensors.direction(id).position];
    sensor_info.t_dir = [sensor_info.t_dir atan2(directions(2,id),directions(1,id))*ones(1,nsensors)];
    data = [data; u_meas.direction(id).field];
    
end
