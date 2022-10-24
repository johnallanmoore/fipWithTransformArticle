# fipWithTransformArticle
Code for IJF article A Study of the Effects of Martensite Phase Transformation on Fatigue Indicator Parameters Determined by a Crystal Plasticity Model

# for polygranular mesh 

abaqus job=rndCpTransElem64kCycle cpus=48 user=umatCpPtNoCommon.f
abaqus python post_FIP_CP.py 
python nonlocalFIP.py  

# for inclusion mesh
abaqus job=fullInclusionXCPTransCycle cpus=48 user=umatCpPtNoCommon.f
abaqus python post_FIP_CP.py 
python nonlocalFIP.py  

# change input names in post_FIP_CP.py and nonlocalFIP.py accordingly. 
# umatCpPtNoCommon.f is avalable on request from anderson.1@osu.edu
# or eamil john.a.moore@marquette.edu