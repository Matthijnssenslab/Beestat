import pandas as pd
castdic = {}
with open('Virdf_metadata.csv') as f:
    for line in f:
        if line.strip().split(',')[1] == 'Sample':
            header = line.strip()
        else:
            sample = line.strip().split(',')[1]
            virus = line.strip().split(',')[2]
            load = line.strip().split(',')[3]
            if sample not in castdic:
                status = line.strip().split(',')[5]
                race = line.strip().split(',')[6]
                Mites =line.strip().split(',')[10]
                txpre =line.strip().split(',')[11]
                txwin =line.strip().split(',')[12]
                txpost =line.strip().split(',')[13]
                castdic[sample] = {'status':status,'mites':Mites}
                if race == 'Carnica':
                    castdic[sample]['Carnica'] = 1
                    castdic[sample]['Unknown'] = 0
                    castdic[sample]['Buckfast + Carnica'] = 0
                    castdic[sample]['Hybrid'] = 0
                    castdic[sample]['Buckfast'] = 0
                    castdic[sample]['Swarm'] = 0
                    castdic[sample]['Queens'] = 0
                elif race == 'Unknown':
                    castdic[sample]['Carnica'] = 0
                    castdic[sample]['Unknown'] = 1
                    castdic[sample]['Buckfast + Carnica'] = 0
                    castdic[sample]['Hybrid'] = 0
                    castdic[sample]['Buckfast'] = 0
                    castdic[sample]['Swarm'] = 0
                    castdic[sample]['Queens'] = 0
                elif race == 'Buckfast + Carnica':
                    castdic[sample]['Carnica'] = 0
                    castdic[sample]['Unknown'] = 0
                    castdic[sample]['Buckfast + Carnica'] = 1
                    castdic[sample]['Hybrid'] = 0
                    castdic[sample]['Buckfast'] = 0
                    castdic[sample]['Swarm'] = 0
                    castdic[sample]['Queens'] = 0
                elif race == 'Hybrid':
                    castdic[sample]['Carnica'] = 0
                    castdic[sample]['Unknown'] = 0
                    castdic[sample]['Buckfast + Carnica'] = 0
                    castdic[sample]['Hybrid'] = 1
                    castdic[sample]['Buckfast'] = 0
                    castdic[sample]['Swarm'] = 0
                    castdic[sample]['Queens'] = 0 
                elif race == 'Buckfast':
                    castdic[sample]['Carnica'] = 0
                    castdic[sample]['Unknown'] = 0
                    castdic[sample]['Buckfast + Carnica'] = 0
                    castdic[sample]['Hybrid'] = 0
                    castdic[sample]['Buckfast'] = 1
                    castdic[sample]['Swarm'] = 0
                    castdic[sample]['Queens'] = 0
                elif race == 'Swarm':
                    castdic[sample]['Carnica'] = 0
                    castdic[sample]['Unknown'] = 0
                    castdic[sample]['Buckfast + Carnica'] = 0
                    castdic[sample]['Hybrid'] = 0
                    castdic[sample]['Buckfast'] = 0
                    castdic[sample]['Swarm'] = 1
                    castdic[sample]['Queens'] = 0
                elif race == 'Queens':
                    castdic[sample]['Carnica'] = 0
                    castdic[sample]['Unknown'] = 0
                    castdic[sample]['Buckfast + Carnica'] = 0
                    castdic[sample]['Hybrid'] = 0
                    castdic[sample]['Buckfast'] = 0
                    castdic[sample]['Swarm'] = 0
                    castdic[sample]['Queens'] = 1
                else:
                    print('problem with race: '+str(race))
                if txpre == 'None':
                    castdic[sample]['PRE_None'] = 1
                    castdic[sample]['PRE_Thymol'] = 0
                    castdic[sample]['PRE_OA'] = 0
                    castdic[sample]['PRE_Fluvalinate'] =0
                    castdic[sample]['PRE_MO'] = 0
                    castdic[sample]['PRE_Coumaphos'] = 0
                    castdic[sample]['PRE_Amitraz'] = 0
                elif txpre == 'Thymol':
                    castdic[sample]['PRE_None'] = 0
                    castdic[sample]['PRE_Thymol'] = 1
                    castdic[sample]['PRE_OA'] = 0
                    castdic[sample]['PRE_Fluvalinate'] =0
                    castdic[sample]['PRE_MO'] = 0
                    castdic[sample]['PRE_Coumaphos'] = 0
                    castdic[sample]['PRE_Amitraz'] = 0
                elif txpre == 'OA':
                    castdic[sample]['PRE_None'] = 0
                    castdic[sample]['PRE_Thymol'] = 0
                    castdic[sample]['PRE_OA'] = 1
                    castdic[sample]['PRE_Fluvalinate'] =0
                    castdic[sample]['PRE_MO'] = 0
                    castdic[sample]['PRE_Coumaphos'] = 0
                    castdic[sample]['PRE_Amitraz'] = 0
                elif txpre == 'Thymol+Fluvalinate':
                    castdic[sample]['PRE_None'] = 0
                    castdic[sample]['PRE_Thymol'] = 1
                    castdic[sample]['PRE_OA'] = 0
                    castdic[sample]['PRE_Fluvalinate'] =1
                    castdic[sample]['PRE_MO'] = 0
                    castdic[sample]['PRE_Coumaphos'] = 0
                    castdic[sample]['PRE_Amitraz'] = 0
                elif txpre == 'MO':
                    castdic[sample]['PRE_None'] = 0
                    castdic[sample]['PRE_Thymol'] = 0
                    castdic[sample]['PRE_OA'] = 0
                    castdic[sample]['PRE_Fluvalinate'] =0
                    castdic[sample]['PRE_MO'] = 1
                    castdic[sample]['PRE_Coumaphos'] = 0
                    castdic[sample]['PRE_Amitraz'] = 0
                elif txpre == 'OA+Thymol':
                    castdic[sample]['PRE_None'] = 0
                    castdic[sample]['PRE_Thymol'] = 1
                    castdic[sample]['PRE_OA'] = 1
                    castdic[sample]['PRE_Fluvalinate'] =0
                    castdic[sample]['PRE_MO'] = 0
                    castdic[sample]['PRE_Coumaphos'] = 0
                    castdic[sample]['PRE_Amitraz'] = 0
                elif txpre == 'OA+Amitraz':
                    castdic[sample]['PRE_None'] = 0
                    castdic[sample]['PRE_Thymol'] = 0
                    castdic[sample]['PRE_OA'] = 1
                    castdic[sample]['PRE_Fluvalinate'] =0
                    castdic[sample]['PRE_MO'] = 0
                    castdic[sample]['PRE_Coumaphos'] = 0
                    castdic[sample]['PRE_Amitraz'] = 1
                elif txpre == 'Thymol+MO':
                    castdic[sample]['PRE_None'] = 0
                    castdic[sample]['PRE_Thymol'] = 1
                    castdic[sample]['PRE_OA'] = 0
                    castdic[sample]['PRE_Fluvalinate'] =0
                    castdic[sample]['PRE_MO'] = 1
                    castdic[sample]['PRE_Coumaphos'] = 0
                    castdic[sample]['PRE_Amitraz'] = 0
                elif txpre == 'OA+Coumaphos':
                    castdic[sample]['PRE_None'] = 0
                    castdic[sample]['PRE_Thymol'] = 0
                    castdic[sample]['PRE_OA'] = 1
                    castdic[sample]['PRE_Fluvalinate'] =0
                    castdic[sample]['PRE_MO'] = 0
                    castdic[sample]['PRE_Coumaphos'] = 1
                    castdic[sample]['PRE_Amitraz'] = 0
                elif txpre == 'Thymol+Coumaphos':
                    castdic[sample]['PRE_None'] = 0
                    castdic[sample]['PRE_Thymol'] = 1
                    castdic[sample]['PRE_OA'] = 0
                    castdic[sample]['PRE_Fluvalinate'] =0
                    castdic[sample]['PRE_MO'] = 0
                    castdic[sample]['PRE_Coumaphos'] = 1
                    castdic[sample]['PRE_Amitraz'] = 0
                elif txpre == 'OA+Fluvalinate':
                    castdic[sample]['PRE_None'] = 0
                    castdic[sample]['PRE_Thymol'] = 0
                    castdic[sample]['PRE_OA'] = 1
                    castdic[sample]['PRE_Fluvalinate'] =1
                    castdic[sample]['PRE_MO'] = 0
                    castdic[sample]['PRE_Coumaphos'] = 0
                    castdic[sample]['PRE_Amitraz'] = 0
                elif txpre == 'Amitraz':
                    castdic[sample]['PRE_None'] = 0
                    castdic[sample]['PRE_Thymol'] = 0
                    castdic[sample]['PRE_OA'] = 0
                    castdic[sample]['PRE_Fluvalinate'] =0
                    castdic[sample]['PRE_MO'] = 0
                    castdic[sample]['PRE_Coumaphos'] = 0
                    castdic[sample]['PRE_Amitraz'] = 1
                else:
                    print('problem with txpre: '+str(txpre))
                if txwin == 'None':
                    castdic[sample]['WIN_None'] = 1
                    castdic[sample]['WIN_Thymol'] = 0
                    castdic[sample]['WIN_OA'] = 0
                    castdic[sample]['WIN_Heat'] =0
                    castdic[sample]['WIN_MO'] = 0
                    castdic[sample]['WIN_Coumaphos'] = 0
                    castdic[sample]['WIN_Chlorophenyl'] = 0
                elif txwin == 'Thymol':
                    castdic[sample]['WIN_None'] = 0
                    castdic[sample]['WIN_Thymol'] = 1
                    castdic[sample]['WIN_OA'] = 0
                    castdic[sample]['WIN_Heat'] =0
                    castdic[sample]['WIN_MO'] = 0
                    castdic[sample]['WIN_Coumaphos'] = 0
                    castdic[sample]['WIN_Chlorophenyl'] = 0
                elif txwin == 'OA':
                    castdic[sample]['WIN_None'] = 0
                    castdic[sample]['WIN_Thymol'] = 0
                    castdic[sample]['WIN_OA'] = 1
                    castdic[sample]['WIN_Heat'] =0
                    castdic[sample]['WIN_MO'] = 0
                    castdic[sample]['WIN_Coumaphos'] = 0
                    castdic[sample]['WIN_Chlorophenyl'] = 0
                elif txwin == 'Heat':
                    castdic[sample]['WIN_None'] = 0
                    castdic[sample]['WIN_Thymol'] = 0
                    castdic[sample]['WIN_OA'] = 0
                    castdic[sample]['WIN_Heat'] =1
                    castdic[sample]['WIN_MO'] = 0
                    castdic[sample]['WIN_Coumaphos'] = 0
                    castdic[sample]['WIN_Chlorophenyl'] = 0
                elif txwin == 'MO':
                    castdic[sample]['WIN_None'] = 0
                    castdic[sample]['WIN_Thymol'] = 0
                    castdic[sample]['WIN_OA'] = 0
                    castdic[sample]['WIN_Heat'] =0
                    castdic[sample]['WIN_MO'] = 1
                    castdic[sample]['WIN_Coumaphos'] = 0
                    castdic[sample]['WIN_Chlorophenyl'] = 0
                elif txwin == 'Coumaphos':
                    castdic[sample]['WIN_None'] = 0
                    castdic[sample]['WIN_Thymol'] = 0
                    castdic[sample]['WIN_OA'] = 0
                    castdic[sample]['WIN_Heat'] =0
                    castdic[sample]['WIN_MO'] = 0
                    castdic[sample]['WIN_Coumaphos'] = 1
                    castdic[sample]['WIN_Chlorophenyl'] = 0
                elif txwin == 'Clorophenyl':
                    castdic[sample]['WIN_None'] = 0
                    castdic[sample]['WIN_Thymol'] = 0
                    castdic[sample]['WIN_OA'] = 0
                    castdic[sample]['WIN_Heat'] =0
                    castdic[sample]['WIN_MO'] = 0
                    castdic[sample]['WIN_Coumaphos'] = 0
                    castdic[sample]['WIN_Chlorophenyl'] = 1
                elif txwin == 'OA+Thymol':
                    castdic[sample]['WIN_None'] = 0
                    castdic[sample]['WIN_Thymol'] = 1
                    castdic[sample]['WIN_OA'] = 1
                    castdic[sample]['WIN_Heat'] =0
                    castdic[sample]['WIN_MO'] = 0
                    castdic[sample]['WIN_Coumaphos'] = 0
                    castdic[sample]['WIN_Chlorophenyl'] = 0  
                else:
                    print('problem with txwin: '+str(txwin))
                if txpost == 'None':
                    castdic[sample]['POST_None'] = 1
                    castdic[sample]['POST_Thymol'] = 0
                    castdic[sample]['POST_OA'] = 0
                    castdic[sample]['POST_Fluvalinate'] =0
                elif txpost == 'Thymol':
                    castdic[sample]['POST_None'] = 0
                    castdic[sample]['POST_Thymol'] = 1
                    castdic[sample]['POST_OA'] = 0
                    castdic[sample]['POST_Fluvalinate'] =0
                elif txpost == 'OA':
                    castdic[sample]['POST_None'] = 0
                    castdic[sample]['POST_Thymol'] = 0
                    castdic[sample]['POST_OA'] = 1
                    castdic[sample]['POST_Fluvalinate'] =0
                elif txpost == 'Fluvalinate':
                    castdic[sample]['POST_None'] = 0
                    castdic[sample]['POST_Thymol'] = 0
                    castdic[sample]['POST_OA'] = 0
                    castdic[sample]['POST_Fluvalinate'] =1
                elif txpost == 'OA+Thymol':
                    castdic[sample]['POST_None'] = 0
                    castdic[sample]['POST_Thymol'] = 1
                    castdic[sample]['POST_OA'] = 1
                    castdic[sample]['POST_Fluvalinate'] =0
                elif txpost == 'Thymol+Fluvalinate':
                    castdic[sample]['POST_None'] = 0
                    castdic[sample]['POST_Thymol'] = 1
                    castdic[sample]['POST_OA'] = 0
                    castdic[sample]['POST_Fluvalinate'] =1
                else:
                    print('problem with txpost: '+str(txpost))                                                       
            if virus == 'AMFV':
                castdic[sample]['AMFVload'] = load
            elif virus == 'bmlv':
                castdic[sample]['BMLVload'] = load
            elif virus == 'thog':
                castdic[sample]['thogload'] = load
            elif virus == 'thika':
                castdic[sample]['thikaload'] = load
            elif virus == 'Unc':
                castdic[sample]['uncload'] = load
            elif virus == 'vdv':
                castdic[sample]['vdvload'] = load
            else:
                castdic[sample]['dwvload'] = load
cumcounter = {}
for i in castdic:
    count = 0
    for j in ['AMFVload','BMLVload','thogload','thikaload','uncload','vdvload','dwvload']:
        if float(castdic[i][j]) > 1000:
            count += 1
    castdic[i]['Numvir'] = count
test = pd.DataFrame(castdic)
test = test.T
test.to_csv("Virdf_metadata_cast.csv",sep=",")

