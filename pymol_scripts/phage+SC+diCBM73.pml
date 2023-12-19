bg_color white

load phage+1000SC.cif.gz, phage

hide cartoon
show spheres
### hides distant parts of phage to simplify rendering
#show spheres, (bychain phage within 250.0 of chain B)
hide spheres, chain SC* and res 111-122;  # hide the SpyTags

color grey, chain P*
color white, chain SC*;             # colour the SpyCatchers by chain

### # highlight NHS-iodoacetate binding targets
### show sphere, A/49/OG
### show sphere, P177/1/N

# define selections of matched SpyTag residues
load selected_prediction.pdb.gz, diCBM1
select diCBM26, diCBM1 and resid 26 and (backbone and not elem H or name CB)
select spytag122, chain SC51 and resid 122 and (backbone or name CB)
pair_fit diCBM26, spytag122

# define selections of matched SpyTag residues
load selected_prediction.pdb.gz, diCBM2
select diCBM26, diCBM2 and resid 26 and (backbone and not elem H or name CB)
select spytag122, chain SC39 and resid 122 and (backbone or name CB)
pair_fit diCBM26, spytag122


# define selections of matched SpyTag residues
load selected_prediction.pdb.gz, diCBM3
select diCBM26, diCBM3 and resid 26 and (backbone and not elem H or name CB)
select spytag122, chain SC45 and resid 122 and (backbone or name CB)
pair_fit diCBM26, spytag122


# define selections of matched SpyTag residues
load selected_prediction.pdb.gz, diCBM4
select diCBM26, diCBM4 and resid 26 and (backbone and not elem H or name CB)
select spytag122, chain SC12 and resid 122 and (backbone or name CB)
pair_fit diCBM26, spytag122

hide cartoon
show spheres, diCBM* and not res 1-25
select SpyTag, res 111-121 and (chain SC12 or chain SC39 or chain SC45 or chain SC51)
show spheres, SpyTag
color orange, SpyTag or diCBM*
# res 122 from diCBM atoms


rotate y, 90
zoom phage, -125


ray 1600, 1200
png phage+1000SC+77diCBM.png
