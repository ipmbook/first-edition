require(grDevices); 


## set path to ghostscript, and put it in the OS environment so embedFonts can find it 
gpathWin = "c:/Program Files (x86)/gs/gs9.07/bin/gswin32c.exe"; # modify as needed
gpathMac = "/usr/local/bin/gs"; # modify as needed
gpath=ifelse(.Platform$OS.type=="windows",gpathWin,gpathMac); 
Sys.setenv(R_GSCMD = gpath)

doEmbed=function(file) {
   embedFonts(file=file, outfile = paste0(file,".embed.eps"), 
   options="-dEmbedAllFonts=true -dEPSCrop")
}   
