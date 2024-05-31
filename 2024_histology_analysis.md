# Quantifying sperm in histological sections

We are interested in quantifying sperm in histological sections. This is difficult to do programatically because (1) sperm stains at a similar intensity as nuclei and (2) sometimes images do not encompass the entire image, so standardizing by the total number of pixels would be inappropriate. As a work around, I've used free software (GIMP) to label sperm with a new color on a new layer. Then I manually erase nuclei and crop the image so no non-testis tissue is in it.  With GIMP I can then easily get the number of pixels that are sperm and also the total number in the image. 

As an example, sample 1841 is a wildtype XL male. Using a 20X image, 342765 out of 6813045 pixels (0.05 = 5%) were sperm. 

Directions:
```
Open the file in GIMP
In the main menu click Windows > Dockable Dialogs > Histogram
In the main menu click Select > By Color 
Adjust the threshold to ~30 on left hand side
In the main menu click Edit > Copy
In the main menu click Edit > Paste as > New Layer
Click on the New Layer on the right side
On the left, change the foreground color to green by clicking on the upper left square that overlaps another square
In the main menu click Edit > Fill with FG color
In the main menu click Select > None
On the left, click on the eraser tool; change the size to ~50  and the hardness/brush setting to 1, change opacity to 100, change aspect ratio to 0
While you are on the New layer with green, erase any green that is not over sperm
In the main menu click Windows > Dockable Dialogs > Histogram
In the main menu click Select > By Color; the number of pixels for the selected color is in the histogram
Now make a new layer:
In the main menu click Layer > New Layer
Change the foreground color to bright red
select the new layer on the right side
Use the paintbucket tool to fill the new layer with red
In the main menu click Select > By Color; the number of pixels for the selected color is in the histogram
```

This should not be so bad because we only have ~5 images per treatment (even less for XT). But it is important because it will allow us to effectively quantify and compare sperm abundance.
