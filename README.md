# Eyeblink analysis workflow

1 - Separate big csv file with all mice into separate individual files (one mouse = one file).  
2 - Make sure that the columns are in the right order, according to your import script. Example below:  
![Screen Shot 2022-06-28 at 09 58 07](https://user-images.githubusercontent.com/57630300/176125681-ad147e3e-bd64-496c-9c48-d7fb3041bc33.png)  
![Screen Shot 2022-06-28 at 09 55 41](https://user-images.githubusercontent.com/57630300/176125196-965f036d-7961-452a-8cb8-9e3fdc92b531.png)  
3 - Run the import script for every mouse. This will generate a matlab variable group per mouse. This will then be imported in the analysis script.  
4 - Run the analysis script. It can be done per mouse or per group. Make sure to change the time windows definitions depending on the eyeblink protocol you are running.  


