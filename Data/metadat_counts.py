import pandas as pd
cumulcount = {}
with open('Virdf_metadata.csv') as f:
    for line in f:
        if line.strip().split(',')[1] != 'Sample':
            if line.strip().split(',')[1] not in cumulcount:
                sample = line.strip().split(',')[1]
                viralload = line.strip().split(',')[3]
                province = line.strip().split(',')[4]
                status = line.strip().split(',')[5]
                race = line.strip().split(',')[6]
                lat = line.strip().split(',')[7]
                lon = line.strip().split(',')[8]
                commune = line.strip().split(',')[9]
                Mitesper100 = line.strip().split(',')[10]
                txpre = line.strip().split(',')[11]
                txwin = line.strip().split(',')[12]
                txpost = line.strip().split(',')[13]
                year = line.strip().split(',')[14]
                cumulcount[sample] = {'Sample':sample,'Virus':'Cumulative','ViralLoad':float(viralload),'Province':province,'Status':status,'Race':race,'LAT':lat,'LON':lon,'Commune':commune,'Mitesper100':Mitesper100,'TXPRE':txpre,'TXWIN':txwin,'TXPOST':txpost,'Year':year}
            else:
                cumulcount[sample]['ViralLoad'] += float(viralload)
test = pd.DataFrame(cumulcount)
test = test.T
test.to_csv("Virdf_metadata_cumulcount.csv",sep=",")
