import os,sys,re,fileinput,string,shutil

##             Dataset        Name   
datasets = [#["SingleMuonGun_NoPU_ieta27_Oct20", "SingleMuonGun_NoPU_ieta27_Oct20"]
            #["SingleMuonGun_wPU_ieta27_Oct20","SingleMuonGun_wPU_ieta27_Oct20"]
            #["QCD_2018D_NoPU_ieta27_Oct20_1","QCD_2018D_NoPU_ieta27_Oct20_1"]
            #["QCD_2018D_NoPU_ieta27_Oct20_2","QCD_2018D_NoPU_ieta27_Oct20_2"]
            ["QCD_2018D_wPU_ieta27_Oct20","QCD_2018D_wPU_ieta27_Oct20"]
            #["WW","WW"]
            #["WZ","WZ"]
            #["ZZ","ZZ"]
            #["ElData1","ElData1"]    
            #["ElData2","ElData2"]
            #["ElData3","ElData3"]
            #["ElData4","ElData4"]
            #["ElData5","ElData5"]
            #["ElData6","ElData6"]
            #["ElData7","ElData7"]
            #["ElData8","ElData8"]    
            #["ElData9","ElData9"]
            #["ElData10","ElData10"]
            #["ElData11","ElData11"]
            #["ElData12","ElData12"]
            #["ElData13","ElData13"]
            #["ElData14","ElData14"]
            #["ElData15","ElData15"]
            #["ElData16","ElData16"]
            #["2017D_SingleMuon","2017D_SingleMuon"]
    #["2017E_SingleMuon","2017E_SingleMuon"]
            #["2017F_SingleMuon","2017F_SingleMuon"]
            #["2017B_SingleMuon","2017B_SingleMuon"]
            #["2017C_SingleMuon","2017C_SingleMuon"]
            #["MuData6","MuData6"]
            #["MuData7","MuData7"]
            #["MuData8","MuData8"]
            #["MuData9","MuData9"]
            #["MuData10","MuData10"]
            #["MuData11","MuData11"]
            #["MuData12","MuData12"]
            #["MuData13","MuData13"]
            #["MuData14","MuData14"]
            #["MuData15","MuData15"]
            #["MuData16","MuData16"]
]

NSections = 10
readFiles = ""

for data in datasets:
    jobidx = 0
    if ( data[0]=="SingleMuonGun_NoPU_ieta27_Oct20"):
        dataname  = "SingleMuonGun_NoPU_ieta27_Oct20"
        inputfname = "SingleMuonGun_NoPU_ieta27_Oct20.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 5

    
    elif ( data[0]=="SingleMuonGun_wPU_ieta27_Oct20"):
        dataname = "SingleMuonGun_wPU_ieta27_Oct20"
        inputfname = "SingleMuonGun_wPU_ieta27_Oct20.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 5    

    elif ( data[0]=="QCD_2018D_NoPU_ieta27_Oct20_1"):
        dataname = "QCD_2018D_NoPU_ieta27_Oct20_1"
        inputfname = "QCD_2018D_NoPU_ieta27_Oct20_1.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 5

    elif ( data[0]=="QCD_2018D_NoPU_ieta27_Oct20_2"):
        dataname = "QCD_2018D_NoPU_ieta27_Oct20_2"
        inputfname = "QCD_2018D_NoPU_ieta27_Oct20_2.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 5

    elif ( data[0]=="QCD_2018D_wPU_ieta27_Oct20"):
        dataname = "QCD_2018D_wPU_ieta27_Oct20"
        inputfname = "QCD_2018D_wPU_ieta27_Oct20.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WW"):
        dataname = "WW"
        inputfname = "WW.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WZ"):
        dataname = "WZ"
        inputfname = "WZ.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    elif ( data[0]=="ZZ"):
        dataname = "ZZ"
        inputfname = "ZZ.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    
    elif ( data[0]=="ElData1"):
        dataname = "ElData1"
        inputfname = "Data1_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="ElData2"):
        dataname = "ElData2"
        inputfname = "Data2_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="ElData3"):
        dataname = "ElData3"
        inputfname = "Data3_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData4"):
        dataname = "ElData4"
        inputfname = "Data4_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData5"):
        dataname = "ElData5"
        inputfname = "Data5_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData6"):
        dataname = "ElData6"
        inputfname = "Data6_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData7"):
        dataname = "ElData7"
        inputfname = "Data7_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ElData8"):
        dataname = "ElData8"
        inputfname = "Data8_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="ElData9"):
        dataname = "ElData9"
        inputfname = "Data9_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="ElData10"):
        dataname = "ElData10"
        inputfname = "Data10_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData11"):
        dataname = "ElData11"
        inputfname = "Data11_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData12"):
        dataname = "ElData12"
        inputfname = "Data12_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData13"):
        dataname = "ElData13"
        inputfname = "Data13_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData14"):
        dataname = "ElData14"
        inputfname = "Data14_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    NFilesTotal = len(readFiles)
    TotalFiles = NFilesTotal
    print "Dataset ",  data[0], " NFilesTotal ", NFilesTotal
    NFilesDone  = 0

    outDir="/eos/cms/store/group/dpg_hcal/comm_hcal/idutta/ntuples/20Oct2018/HCALSample"+data[0]
    print outDir
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    while( NFilesDone < NFilesTotal ) :
        thisList = readFiles[NFilesDone : NFilesDone+NSections]
        print "NFilesDone ", NFilesDone, "len(thisList)", len(thisList)

        ##you may have to give full path i.e. CurrentDIR/condor_submit/runlist_...
        inputRunListName = "/afs/cern.ch/work/i/idutta/private/HCAL_muALCA/CMSSW_10_1_6/src/Muon_Calibration_noTrack_2018/muon_calibration_ana/bash/condor/condor_submit/runList_"+data[0]+"_"+str(jobidx)+".txt"
        inputRunList = open(inputRunListName, "w")
        for line in thisList:
            inputRunList.write(line)

        condorSubmit = "condor_submit/submitCondor_"+data[0]+"_"+str(jobidx)
        jobName      = "20Oct2018"+data[0]+"_job"+str(jobidx)
        #outHistFile = data[0]+"_job"+str(jobidx)+".root"
        #isData       ="T"
        #isData       ="F"
        shutil.copyfile("proto_condor_submit",condorSubmit)
        for line in fileinput.FileInput(condorSubmit, inplace=1):
            line=line.replace("JOBNAME", jobName)
            line=line.replace("FILELIST",inputRunListName)
            #line=line.replace("ROOTOUT",outHistFile)
            #line=line.replace("DATANAME",dataname)
            #line=line.replace("ISDATA",isData)
            line=line.replace("OUTDIR",outDir)
            print line.rstrip()
        
        submitCommand = "condor_submit "+condorSubmit
        print submitCommand
        os.system(submitCommand)     
        jobidx = jobidx+1
        NFilesDone = NFilesDone + len(thisList)

    print "Final NFilesDone ", NFilesDone
