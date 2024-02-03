# -*- coding: utf-8 -*-
#%%
"""
Created on Thu May 11 10:51:00 2023

@author: ganderson
"""

import geopandas
import yaml
import numpy as np
import logging
from shapely.geometry import Point
import os
import pandas as pd

try:
    os.remove('NetworkGen.log')
except:
    pass


logging.basicConfig(filename='NetworkGen.log', level=logging.WARNING, filemode='w',
                    format='%(asctime)s:%(levelname)s:%(lineno)d:%(message)s')


# File = 'Shapefiles/BessemerPointLayerForConversion.shp'
File = 'D:/project/altrios-public/docs/src/doc/network-files/Shapefiles/BessemerTry2PointData.shp'

# LineFile = 'Shapefiles/BessemerLineLayerForConversion.shp'
LineFile ='D:/project/altrios-public/docs/src/doc/network-files/Shapefiles/BessemerTry2LineData.shp'

SpeedRestrionFile = 'Other Data/BessemerSpeedRestrictions.csv'


# NetworkName = 'BessemerAndLakeErie'
NetworkName = 'BessemerAndLakeErieForVideo2'


logging.debug('Network: '.format(NetworkName))
logging.debug('Point File: '.format(File))
logging.debug('Line File: '.format(LineFile))

LimitType = {1 : 'CivilSpeed', 
             2 : 'MaxPermissibleSpeed',
             3 : 'MassTotal',
             4 : 'MassPerBrake',
             5 : 'AxleCount',
             6 : 'Unknown Limit Type'}

CompareType = {1 : 'TpEqualRp',
               2 : 'TpGreaterThanRp',
               3 : 'TpLessThanRp',
               4 : 'TpGreaterThanEqualRp',
               5 : 'TpLessThanEqualRp'
               }

TrainType = {1 : 'Freight',
             2 : 'Passenger',
             3 : 'Intermodal',
             4 : 'HighSpeedPassenger',
             5 : 'TiltTrain',
             6 : 'Commuter'}


# %%

def BlockToDict(BlockData):
    BlockData = BlockData.sort_values('distance')
    BlockDict = {'idx_next' : 0,
    'idx_next_alt' : 0,
    'idx_prev' : 0,
    'idx_prev_alt' : 0,
    'idx_curr' : 0,
    'idx_flip' : 0,
    'osm_id' :  BlockData.osm_id.values[0],
    # 'LinkedToNetwork' : False, #dont think i need because its being tracked in the data dataframe.
    # 'subidv' : BlockData.SUBDIV.values[0]
    }
    

    Elevs = []
    Headings = []
    for idx, row in BlockData.iterrows():
        Elevs.append({'offset': row['distance'],
                      'elev' : row['Elevation'],
                      })
        
        Headings.append({'offset': row['distance'],
                      'heading' : row['angle']*2*np.pi/360,
                      'Lat' : list(row.geometry.coords)[0][0],
                      'Lon' : list(row.geometry.coords)[0][0],
                      'Elev' : row['Elevation']})
        
    BlockDict['elevs'] = Elevs
    BlockDict['headings'] = Headings
    BlockDict['length'] = float(BlockData['distance'].max())
    # try:
        
    SpeedSets = []    
        
    Subdiv = BlockData.name.mode()
    
    if Subdiv.shape[0] > 0:
        SpeedRestriction = SpeedRestrictions[SpeedRestrictions.Subdivision == Subdiv.values[0]]
    else:
        SpeedRestriction = pd.DataFrame()
    
    if SpeedRestriction.shape[0] > 0:
        MaxSpeed = float(SpeedRestriction.Limit.values[0] / 2.23693629) #TODO: set limits through sidings and turnoouts, limits based on milepost, limits based on TPOB
    else:
        MaxSpeed = 6.7056 #15 mph
    
    for key in TrainType.keys():
          Speeds = [{'offset_start': 0,
                    'offset_end' : float(BlockData['distance'].max()),
                    'speed' : MaxSpeed}]
         
          SpeedDict = {'speed_limits': Speeds,
                        'speed_params': [], #{'limit_val': float(row1.RestrictionParameter),
                      #                  'limit_type': LimitType[row1.RestrictionType],
                      #                  'compare_type': CompareType[row1.RestrictionParamOperator]},
                      # can fix this if we omit line 168 and add some more smarts
                      'train_type': TrainType[key],
                      'is_head_end': False}
         
          SpeedSets.append(SpeedDict)

    BlockDict['speed_sets'] = SpeedSets
    
    if BlockData.osm_id.values[0].find('_reverse') >-1:  
        # if true this is a reverse block;  the next link will be contained in the idx_next_reverse, and the previous links will be in the idx_next.  Forward links are the opposite of this.
        PrevSuffix = ''
        NextSuffix = '_reverse'
    else:
        PrevSuffix = '_reverse'
        NextSuffix = ''
    
    BlockDict['idx_next'] = BlockData['idx_next' + NextSuffix].values[0]
    BlockDict['idx_next_alt'] = BlockData['idx_next_alt' + NextSuffix].values[0]
    BlockDict['idx_prev'] = GetFlippedOSMID(BlockData['idx_next' + PrevSuffix].values[0])
    BlockDict['idx_prev_alt'] = GetFlippedOSMID(BlockData['idx_next_alt' + PrevSuffix].values[0])
        

    return BlockDict

def BaseOSMExtract(OSMString):
    #need this function to take the _reverse off the osm_id to select the right line from the LineData that doesn't need to be reversed.
    return OSMString.replace('_reverse','')

def GetFlippedOSMID(OSMString):#should probably rename to something other than string.
    if OSMString == 0:
        return 0
    if OSMString == BaseOSMExtract(OSMString):
        return OSMString + '_reverse' #TODO:  turn this into a one line function so _reverse is only in one spot or atleast make a global variable.
    else:
        return BaseOSMExtract(OSMString)
        

def LinkUpBlock(BlockDict, BlockData, LineData):
    
    StartRow = BlockData[BlockData['distance']==0].copy() #grab end of line where is starts
    StartRow['base_osm_id'] = StartRow.osm_id.apply(BaseOSMExtract)
    StartBuff = StartRow.buffer(1.0) #create 2m buffer off of point.  Needs to be in correct crs
    StartBuff =geopandas.GeoDataFrame(geometry=StartBuff) #convert series to dataframe
    StartLinks = geopandas.sjoin( Data, StartBuff) #join complete set of points to buffer to find connecting blocks
    StartLinks['base_osm_id'] = StartLinks.osm_id.apply(BaseOSMExtract)
    StartLinks = StartLinks[StartLinks.base_osm_id != StartRow.base_osm_id.values[0]] #drop the link that we're currently linking up
    
    
    EndRow = BlockData[BlockData['distance'] == BlockData['distance'].max()].copy()
    EndRow['base_osm_id'] = EndRow.osm_id.apply(BaseOSMExtract)
    EndBuff = EndRow.buffer(1.0) #create 2m buffer off of point.  Needs to be in correct crs
    EndBuff =geopandas.GeoDataFrame(geometry=EndBuff) #convert series to dataframe
    EndLinks = geopandas.sjoin( Data, EndBuff) #join complete set of points to buffer to find connecting blocks
    EndLinks['base_osm_id'] = EndLinks.osm_id.apply(BaseOSMExtract)
    EndLinks = EndLinks[EndLinks.base_osm_id != EndRow.base_osm_id.values[0]] #drop the link that we're currently linking up

    #creating buffer from line of this section of track, grabbing all lines that are in link lists, and then finding which ones cross the buffers.  This will eliminate turnouts that are not able to be turned onto
    
    #creating 13 foot buffere based on NS track requirement that tracks be 14 feet apart.  Anything that crosses this buffer can't be a potential link and will be blocked out by this section of track. 
    #http://nscorp.com/content/dam/nscorp/industrial-development/track-design-information/NS-Standards-for-Industry-Tracks.pdf
    BlockBuffer = LineData[LineData.osm_id == StartRow.base_osm_id.values[0]].buffer(3.96,cap_style=2) #grab the line to create the buffer
    BlockBuffer=geopandas.GeoDataFrame(geometry=BlockBuffer)
    LinkLines = LineData[LineData.osm_id.isin(StartLinks.base_osm_id.tolist() + EndLinks.base_osm_id.tolist())]
    BadLinkLines = geopandas.sjoin(LinkLines, BlockBuffer, predicate='crosses') #TODO: use this for blocking out sections of track
    
    StartLinks = StartLinks[~StartLinks.osm_id.isin(BadLinkLines.osm_id)]
    EndLinks = EndLinks[~EndLinks.osm_id.isin(BadLinkLines.osm_id)]
    
    MakeTroubleShootingPlots = False
    if MakeTroubleShootingPlots:
        ax=BlockBuffer.plot(color='red')
        LinkLines.plot(ax=ax, color=['g','b','y','m'])
    
    StartLinks['Heading Error'] = np.mod(StartLinks.angle - StartRow.angle.values[0], 180)
    StartLinks.sort_values('Heading Error', inplace=True)
    
    if StartRow.osm_id.values[0].find('_reverse') > -1:
        Suffix = '_reverse'
    else:
        Suffix = ''
    
    if StartLinks.shape[0] > 0:
        BlockDict['idx_prev'] = StartLinks.osm_id.values[0] + Suffix
        BlockDict['idx_prev_osm_id'] = StartLinks.osm_id.values[0] + Suffix
    if StartLinks.shape[0] > 1:
        BlockDict['idx_prev_alt'] = StartLinks.osm_id.values[1] + Suffix
        BlockDict['idx_prev_alt_osm_id'] = StartLinks.osm_id.values[1] + Suffix
    if StartLinks.shape[0] > 2:
        assert 'Bad things are going down with too many link connecting at start end osm_id: {}'.format(StartRow.osm_id.values[0])

    EndLinks['Heading Error'] = np.mod(EndLinks.angle - EndRow.angle.values[0], 180)
    EndLinks.sort_values('Heading Error', inplace=True)
    
    if EndLinks.shape[0] > 0:
        BlockDict['idx_next'] = EndLinks.osm_id.values[0] + Suffix
        BlockDict['idx_next_osm_id'] = EndLinks.osm_id.values[0] + Suffix
    if EndLinks.shape[0] > 1:
        BlockDict['idx_next_alt'] = EndLinks.osm_id.values[1] + Suffix
        BlockDict['idx_next_alt_osm_id'] = EndLinks.osm_id.values[1] + Suffix
    if EndLinks.shape[0] > 2:
        assert 'Bad things are going down with too many link connecting at End end osm_id: {}'.format(EndRow.osm_id.values[0])
   
    return BlockDict
    
def FlipBlock(BlockData):
    
    FlippedData = BlockData.copy()
    FlippedData.angle = (FlippedData.angle + 180) % 360
    FlippedData['distance'] = ((FlippedData['distance'] - FlippedData['distance'].max()) * -1).abs()
    FlippedData.sort_values('distance', inplace=True)
    FlippedData['osm_id'] = FlippedData['osm_id'] + '_reverse'

    FlippedDict = BlockToDict(FlippedData)
    return FlippedDict, FlippedData



# def BlockReverse(DataSeries):
#     # print(DataSeries.index)
#     ReversedSeries = DataSeries.copy()
#     ReversedSeries.Grade = ReversedSeries.Grade * -1
#     ReversedSeries.Heading = (ReversedSeries.Heading + 180) % 360
#     ReversedSeries.IDNext = DataSeries.IDPrev
#     ReversedSeries.IDNextAlt = DataSeries.IDPrevAlt
#     ReversedSeries.IDPrev = DataSeries.IDNext
#     ReversedSeries.IDPrevAlt = DataSeries.IDNextAlt
    
#     if ReversedSeries.BlockDirection == 'Ascending':
#         ReversedSeries.BlockDirection = 'Descending'
#         ReversedSeries.UID = DataSeries.UID + '_Desc'
#     else:
#         ReversedSeries.BlockDirection = 'Ascending'
#         ReversedSeries.UID = DataSeries.UID.replace('_Desc', '')
    
#     return ReversedSeries   

def FlipOSMID(osm_id):
    if osm_id.find('_reverse') > -1:
        return osm_id.replace('_reverse')
    else:
        return osm_id  + '_reverse'
     


def TraverseNextLink(BlockData, Data, LineData, BlockList):

    return Data, BlockList



#%%

Data = geopandas.read_file(File)
Data = Data.to_crs('ESRI:102009')
Data.osm_id = Data.osm_id_plu
LineData = geopandas.read_file(LineFile)
LineData = LineData.to_crs('ESRI:102009')
LineData.osm_id = LineData.osm_id_plu
UniqueIDs = Data.osm_id.unique()

StillLinksToConnect = True
LinkLengths =  Data.groupby('osm_id')['distance'].max().to_frame()
LinkLengths['osm_id'] = LinkLengths.index.values
LinkLengths = LinkLengths.reindex()
Data = Data.join(LinkLengths, on='osm_id', rsuffix='Max')
Data['idx_next_osm_id'] = 0
Data['idx_next_alt_osm_id'] = 0
Data['idx_prev_osm_id'] = 0
Data['idx_prev_alt_osm_id'] = 0
Data['idx_next_osm_id_reverse'] = 0
Data['idx_next_alt_osm_id_reverse'] = 0
Data['idx_prev_osm_id_reverse'] = 0
Data['idx_prev_alt_osm_id_reverse'] = 0
Data['Potential Links_reverse'] = 0
Data['Mod. Heading'] = np.mod(Data['angle'], 180)

SpeedRestrictions = pd.read_csv(SpeedRestrionFile)

#this loop will go through the code and find all possible links within the data
for UID in UniqueIDs:
    print(UID)
    BlockData = Data[Data.osm_id==UID]
    BlockData = BlockData.sort_values('distance')
    
    StartRow = BlockData[BlockData['distance']==0].copy() #grab end of line where is starts
    StartRow['base_osm_id'] = StartRow.osm_id.apply(BaseOSMExtract)
    StartBuff = StartRow.buffer(.05) #create 2m buffer off of point.  Needs to be in correct crs
    StartBuff =geopandas.GeoDataFrame(geometry=StartBuff) #convert series to dataframe
    StartLinks = geopandas.sjoin( Data, StartBuff) #join complete set of points to buffer to find connecting blocks
    StartLinks['base_osm_id'] = StartLinks.osm_id.apply(BaseOSMExtract)
    StartLinks = StartLinks[StartLinks.base_osm_id != StartRow.base_osm_id.values[0]] #drop the link that we're currently linking up
    
    
    EndRow = BlockData[BlockData['distance'] == BlockData['distance'].max()].copy()
    EndRow['base_osm_id'] = EndRow.osm_id.apply(BaseOSMExtract)
    EndBuff = EndRow.buffer(.05) #create 2m buffer off of point.  Needs to be in correct crs
    EndBuff =geopandas.GeoDataFrame(geometry=EndBuff) #convert series to dataframe
    EndLinks = geopandas.sjoin( Data, EndBuff) #join complete set of points to buffer to find connecting blocks
    EndLinks['base_osm_id'] = EndLinks.osm_id.apply(BaseOSMExtract)
    EndLinks = EndLinks[EndLinks.base_osm_id != EndRow.base_osm_id.values[0]] #drop the link that we're currently linking up

    #creating buffer from line of this section of track, grabbing all lines that are in link lists, and then finding which ones cross the buffers.  This will eliminate turnouts that are not able to be turned onto
    
    #creating 13 foot buffere based on NS track requirement that tracks be 14 feet apart.  Anything that crosses this buffer can't be a potential link and will be blocked out by this section of track. 
    #http://nscorp.com/content/dam/nscorp/industrial-development/track-design-information/NS-Standards-for-Industry-Tracks.pdf
    BlockBuffer = LineData[LineData.osm_id == StartRow.base_osm_id.values[0]].buffer(3.96,cap_style=2) #grab the line to create the buffer
    BlockBuffer=geopandas.GeoDataFrame(geometry=BlockBuffer)
    LinkLines = LineData[LineData.osm_id.isin(StartLinks.base_osm_id.tolist() + EndLinks.base_osm_id.tolist())]
    BadLinkLines = geopandas.sjoin(LinkLines, BlockBuffer, predicate='crosses') #TODO: use this for blocking out sections of track
    BadLinkLines = pd.concat([BadLinkLines, geopandas.sjoin(LinkLines, BlockBuffer, predicate='within')])
    MakeTroubleShootingPlots = False
    if MakeTroubleShootingPlots:
        import folium
        ax=BlockBuffer.plot(color='red')
        LinkLines.plot(ax=ax, color=['g','b','y','m'])
        
        # LinkMap=LinkLines.explore().save('LinkMap.html')
        
    
    StartLinks = StartLinks[~StartLinks.osm_id.isin(BadLinkLines.osm_id)]
    EndLinks = EndLinks[~EndLinks.osm_id.isin(BadLinkLines.osm_id)]
    
    StartLinks['Heading Error'] = (StartLinks['Mod. Heading'] - StartRow['Mod. Heading'].values[0]).abs()
    EndLinks['Heading Error'] = (EndLinks['Mod. Heading'] - EndRow['Mod. Heading'].values[0]).abs()
    StartLinks = StartLinks[StartLinks['Heading Error'] < 30] #this filters out the tracks that cross at right angles
    EndLinks = EndLinks[EndLinks['Heading Error'] < 30] #this filters out the tracks that cross at right angles
    
    if StartLinks.shape[0] > 0 :
        Data.loc[(Data.osm_id == StartRow.osm_id.values[0]) & (Data['distance']==0),['Potential Links']] = str(StartLinks.osm_id.to_list())#remember that the current id will show up in the link that it is matching to get the exact poitn
        Data.loc[(Data.osm_id == StartRow.osm_id.values[0]) & (Data['distance']==0),['Potential Links Index']] = str(StartLinks.index.to_list())#remember that the current id will show up in the link that it is matching to get the exact poitn

    else:
        Data.loc[(Data.osm_id == StartRow.osm_id.values[0]) & (Data['distance']==0),['Potential Links']] = 'None'
        Data.loc[(Data.osm_id == StartRow.osm_id.values[0]) & (Data['distance']==0),['Potential Links Index']] = 'None'
    if EndLinks.shape[0] > 0 :
        Data.loc[(Data.osm_id == EndRow.osm_id.values[0]) & (Data['distance']==Data['distanceMax']),['Potential Links']] = str(EndLinks.osm_id.to_list())
        Data.loc[(Data.osm_id == EndRow.osm_id.values[0]) & (Data['distance']==Data['distanceMax']),['Potential Links Index']] = str(EndLinks.index.to_list())
    
    else:
        Data.loc[(Data.osm_id == EndRow.osm_id.values[0]) & (Data['distance']==Data['distanceMax']),['Potential Links']] = 'None'
        Data.loc[(Data.osm_id == EndRow.osm_id.values[0]) & (Data['distance']==Data['distanceMax']),['Potential Links Index']] = 'None'


#%%

def GenConnectionOSMIDWithSuffix(SeriesData):
    if SeriesData.distance == 0:
        return ''
    else:
        return '_reverse'


point = Point(0, 0)
Data['Distance'] = Data.distance(point)
Data['Traversed'] = False
Data['Traversed_reverse'] = False        
LineEnds = Data[Data['Potential Links']=='None'].copy()
LineEnds = LineEnds.sort_values('Distance', ascending=False)

LinkEnds = Data[Data['Potential Links'].notna()] #selecting ends of links so that each link can be traversed both directions
rng = np.random.default_rng(12345)


# for idx, row in LineEnds.iterrows(): #this was to iterate from only end points.  Commented out to iterate through from each link
for idx, row in LinkEnds.iterrows():    
    print(idx)
    PathComplete = False
    # CurrentBlock = {'osm_id' : row['osm_id']}
    # CurrentBlock['Forward'] = row.distance == 0
    CurrentBlock = row  
    CurrentPathLength = 1
    logging.debug("-------------------------------------------------------------------------------------")
    logging.debug('Start Link: {}, Start Distance: {}'.format(CurrentBlock.osm_id, CurrentBlock.distance))
    PathLinks = []
    while not PathComplete:
        #current link suffix is handling if we're going through the link forward or reverse
        CurrentLinkSuffix = GenConnectionOSMIDWithSuffix(CurrentBlock)
        Data.loc[Data.osm_id == CurrentBlock.osm_id, ['Traversed' + CurrentLinkSuffix]] = True #this will set the proper direction for the link to traversed.
        PathLinks.append(CurrentBlock.osm_id + CurrentLinkSuffix)
        #this grab the 0 offset or the max offset point for the link.  If we start at 0, we're going forward.
        #We need to grab the max offset point that is on the other end to link to the next link.
        #If we are going in reverse we start at the max end and need to grab the 0 end.
        EndPointData = Data[(Data.osm_id == CurrentBlock['osm_id'])].sort_values(['distance'])
        if CurrentLinkSuffix == '':
            EndPointData = EndPointData.iloc[-1,:]
        else:
            EndPointData = EndPointData.iloc[0,:]

        
        #This grabs the points that over
        
        if EndPointData['Potential Links Index'] != 'None':
            PotentialLinks = list(map(int,EndPointData['Potential Links Index'][1:-1].replace("'","").split(','))) #getting all the indexes out as integers here.  Just takes a few steps
        else:
            PotentialLinks = pd.DataFrame()
        # if len(PotentialLinks) > 1: #this is a intended error to catch a multilink link for developing code.  Need to delete or comment out at some point.
        #         x=asdf
        PotentialLinks = Data[Data.index.isin(PotentialLinks)].copy() #this should only be the row that is the end of the link that overlaps with the current link.  Longest this should be is two.
        
        #calculate the difference heading to figure out which is alt and which is normal
        PotentialLinks['Heading Error'] = np.mod(PotentialLinks.angle - EndPointData.angle, 180)
        PotentialLinks.sort_values('Heading Error', inplace=True)
        


        Data['Complete']=Data['Traversed'] & Data['Traversed_reverse']
        NetworkDone = Data['Complete'].min()
        if PotentialLinks.shape[0] > 0:
            #grab the closest heading and set it to the next link
            Suffix1 = GenConnectionOSMIDWithSuffix(PotentialLinks.iloc[0,:])
            Data.loc[(Data.osm_id == CurrentBlock['osm_id']), ['idx_next_osm_id'+ CurrentLinkSuffix]] = PotentialLinks.osm_id.values[0] + Suffix1
            Data.loc[(Data.osm_id == CurrentBlock['osm_id']), ['idx_next'+ CurrentLinkSuffix]] = PotentialLinks.osm_id.values[0] + Suffix1
        if PotentialLinks.shape[0] > 1:
            #grab the other link and set it to next alt.
            Suffix2 = GenConnectionOSMIDWithSuffix(PotentialLinks.iloc[1,:])
            Data.loc[(Data.osm_id == CurrentBlock['osm_id']), ['idx_next_alt_osm_id'+ CurrentLinkSuffix]] = PotentialLinks.osm_id.values[1] + Suffix2
            Data.loc[(Data.osm_id == CurrentBlock['osm_id']), ['idx_next_alt'+ CurrentLinkSuffix]] = PotentialLinks.osm_id.values[1] + Suffix2
        if PotentialLinks.shape[0] > 2:
            #if you have three links, there is a problem with the network..... May the code.
            assert 'Bad things are going down with too many link connecting at start end osm_id: {}'.format(StartRow.osm_id.values[0])
        
        
        if (PotentialLinks.shape[0] == 0) or (CurrentPathLength > 250): #len(UniqueIDs)):
            #if not links, we hit the end and should stop this path
            PathComplete=True
            
            if CurrentPathLength > 2:
                logging.warning('Path length is quite long.  May indicate loop in track that is not correct')
                logging.warning('Link Path: {}'.format(PathLinks))
                
            logging.debug('Complete Path Length: {}'.format(CurrentPathLength))
            logging.debug('******')
           
        else:
            CurrentPathLength += 1
            
            if PotentialLinks['Traversed'+Suffix1].values[0] == False:
                CurrentBlock = PotentialLinks.iloc[0,:]
            elif PotentialLinks['Traversed'+Suffix1].values[0] == False:
                CurrentBlock = PotentialLinks.iloc[1,:]
            else:
                RandomPick =np.random.choice([0,1])
                RandomPick = np.min([RandomPick, PotentialLinks.shape[0]-1])
                CurrentBlock = PotentialLinks.iloc[RandomPick,:]
                
            logging.debug('Next Block: {}, Next Point Offset: {}, Path Length: {}'.format(CurrentBlock.osm_id, CurrentBlock.distance, CurrentPathLength))
        #need to put logic here to look at traverse status of next link
        #go to normal link by default
        #then go to alt if normal is done
        #then randomly select leg if both checked
        
        
        
        Data['Complete']=Data['Traversed'] & Data['Traversed_reverse']
#%%
BlockList = [{'elevs' : [],
           'speed_sets' : [],
           'idx_next' : 0,
           'idx_next_alt' : 0,
           'idx_prev' : 0,
           'idx_prev_alt' : 0,
           'idx_curr' : 0,
           'idx_flip' : 0,
           'length' : 0,
           'osm_id' : 0}]

FlippedBlockList = []


cols = ['idx_next', 'Complete', 'idx_next_reverse','idx_next_alt_reverse', 'idx_next_alt']
for col in cols:
    Data[col] = Data[col].fillna(0)

for UID in UniqueIDs:
    BlockData = Data[Data.osm_id==UID]
    BlockDict = BlockToDict(BlockData)
    FlippedDict, FlippedData = FlipBlock(BlockData)
    # BlockDict = LinkUpBlock(BlockDict, BlockData, LineData)
    # FlippedDict = LinkUpBlock(FlippedDict, FlippedData, LineData)
    BlockList.append(BlockDict)
    FlippedBlockList.append(FlippedDict)
    #TODO check for link leading back to current link (trap where it gets caught between two links)
    #DONE TODO create dict of idx and osm_id on complete blocklist
    #DONE TODO replace idx_next,..... with integers
    #TODO check heading angle at junctions (low priority)
    #TODO add blocking functionality
    print(UID)
    
BlockList = BlockList + FlippedBlockList





#%% replacing osm_id with integer for index

BlockIDDict = {}
for i in range(len(BlockList)):
    BlockIDDict[BlockList[i]['osm_id']] = i
    BlockList[i]['idx_curr'] = i
    print(i)

idx_keys = ['idx_next', 'idx_next_alt', 'idx_prev', 'idx_prev_alt']#, 'idx_flip']
for i in range(len(BlockList)):
    if BlockList[i]['osm_id'] != 0:
        BlockList[i]['idx_flip'] = BlockIDDict[GetFlippedOSMID(BlockList[i]['osm_id'])]
        for key in idx_keys:
                BlockList[i][key] = BlockIDDict[BlockList[i][key]]  


#%%
with open(NetworkName + '.yaml', 'w') as f:
    f.write(yaml.dump(BlockList))

