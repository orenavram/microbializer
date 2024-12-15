Gallery backups.

If for some reason microbializer's gallery was deleted, say, **GammaProteoBacteria dataset is unavailable. copy the GammaProteoBacteria folder to /var/www/html/microbializer (currently it's on powerweb1, but who knows..). so that afterwards you'll have this (remote) path valid: /var/www/html/microbializer/GammaProteoBacteria
E.g, use:
scp -r /Users/Oren/Dropbox/Projects/microbializerPower/case_studies/GammaProteoBacteria bioseq@powerweb1:/var/www/html/microbializer/


Next, ssh to this remote folder
ssh bioseq@powerweb1.tau.ac.il
cd /var/www/html/microbializer/GammaProteoBacteria

and extract M1CR0B1AL1Z3R_GammaProteoBacteria_outputs.zip using (it will take a few minutes):
unzip -q -d M1CR0B1AL1Z3R_GammaProteoBacteria_outputs/ M1CR0B1AL1Z3R_GammaProteoBacteria_outputs.zip

Then, the gallery of the GammaProteoBacteria should be available at:
https://microbializer.tau.ac.il/GammaProteoBacteria/output.html


** If it's another folder name, just find+replace *all* "GammaProteoBacteria" instances in this file to the folder name you need