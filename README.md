
<h1 align="center">Satellite Tracker</h1>

## About The Project

Project built for my Final Degree Project at the Universitat Politecnica de Catalunya (UPC) in the Double Bachelor's Degree in Aerospace Systems Engineering and Network Engineering. Capable of showing the position of multiple satellites, during their full next orbit around the earth in a 2D map. Also capable of shwowing the observed Doppler effect and when and in which direction will the satellites will be visible for the user. Built entirely by Marc Xapelli with help of professor Juan José Olmos Bonafé. The backend was built in Python and the frontend in Javascript using Visual Studio Code, it uses a MongoDB database. It is deployed in Heroku. To go to the website click[here](https://satellite-tracker-eetac.herokuapp.com). 

## User Guide

![alt text](https://github.com/AsterixDecoder/AsterixDecoder/blob/main/AsterixDecoder/images/loadFile.PNG?raw=true)

When opening, upload a file --> Load File to select the .ast file.

After Loading, open a Table to see the information of the respective category.

![alt text](https://github.com/AsterixDecoder/AsterixDecoder/blob/main/AsterixDecoder/images/cat21.PNG?raw=true)

You can click on the cell that say "Click to expand", to see all the information of a certain Data Item. You can search also for packets using the Track Number in CAT10 and Target Identification in CAT21.

![alt text](https://github.com/AsterixDecoder/AsterixDecoder/blob/main/AsterixDecoder/images/search.PNG?raw=true)

In the Map View you can simulate the flights and there are checks to filter by type of emission and buttons to accelerate or decelerate the time. 
![alt text](https://github.com/AsterixDecoder/AsterixDecoder/blob/main/AsterixDecoder/images/map.PNG?raw=true)

You can also search for a determined Target Identification or Track Number, and look at all the position of a determined flight.
![alt text](https://github.com/AsterixDecoder/AsterixDecoder/blob/main/AsterixDecoder/images/viewOld.PNG?raw=true)