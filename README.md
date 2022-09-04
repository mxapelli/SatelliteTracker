
<h1 align="center">Satellite Tracker</h1>

## About The Project

Project built for my Final Degree Project at the Universitat Politecnica de Catalunya (UPC) in the Double Bachelor's Degree in Aerospace Systems Engineering and Network Engineering. Capable of showing the position of multiple satellites, during their full next orbit around the earth in a 2D map. Also capable of shwowing the observed Doppler effect and when and in which direction will the satellites will be visible for the user. Built entirely by Marc Xapelli with help of professor Juan José Olmos Bonafé. The backend was built in Python and the frontend in Javascript using Visual Studio Code, it uses a MongoDB database. It is deployed in Heroku. To go to the website click [here](https://satellite-tracker-eetac.herokuapp.com). 

## User Guide

![alt text](./images/index.PNG?raw=true)

Select a specific satellite or an entire constellation of satellites. You can also allow the use of your real geolocation.

After selecting what you want to see, a new tab will be opened with the trajectory of the satellite or satellites and the visibility area of the user.

![alt text](./images/constellation.PNG?raw=true)

In the constellation view if you select a satellite a new tab will be opened with the specific information of the satellite.

![alt text](./images/constellation2.PNG?raw=true)

In the Satellite View if the satellite is visible, you can see more data by clicking on the button.
![alt text](./images/sat1.PNG?raw=true)

The Analysis View show more data about the visibility of the satellite and other parameters.
![alt text](./images/satdata.PNG?raw=true)

You can also search for a determined satellite name or catalog number.
![alt text](./images/search.PNG?raw=true)