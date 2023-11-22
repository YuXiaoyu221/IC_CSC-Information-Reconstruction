function rgb_img = WriteMultiBandsImages(rgb_img,name,nBits)
[~, ~, bands] = size(rgb_img);
t = Tiff(name,'w');  %creat tiff

tagstruct.ImageLength = size(rgb_img,1); 
tagstruct.ImageWidth = size(rgb_img,2);  

tagstruct.Photometric = 1;

%unit8
tagstruct.BitsPerSample = nBits;
tagstruct.SamplesPerPixel = bands;
tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB'; 
tagstruct.SampleFormat = 1;
t.setTag(tagstruct);

t.write(rgb_img);
t.close;
end

