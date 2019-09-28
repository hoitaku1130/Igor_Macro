#pragma rtGlobals=3	// Use modern global access method and strict wave access.
//#include <WaveSelectorWidget>
//#include<PopupWaveSelector>

Menu "Noise_Filter"
	"Noise_fil"
	"Make_Absorbance"
	"autoFit_Panel"
	"subtract"
End

static Function /S get_datafolder_from_str(str)
	string str
	if(cmpstr(str,"root") == 0)
		return "root:"
	else
		return  "root:" + str
	endif
End

Function make_zero(dfr, str , offset)
	DFREF dfr
	String str
	Variable offset
	Variable i
	Wave /SDFR=dfr my_wave = $(str)
	variable max_of_index = numpnts(my_wave)
	for (i=offset;i<max_of_index;i+=1)
		my_wave[i] = 0
	endfor
End

static Function cutting_data_by_step_func(dfr,str,offset,sigma)
	DFREF dfr
	String str
	variable offset,sigma
	variable i
	wave /C /SDFR=dfr w = $(str)
	variable max_of_index = numpnts(w)
	for(i=0;i<max_of_index;i+=1)
		variable tmp_v = NF_step_func(sigma,offset,i) 
		w[i] = cmplx(real(w[i])*tmp_v,imag(w[i])*tmp_v)
	endfor
End

static Function NF_step_func(sig,off_x,x)
	variable sig,off_x,x
	variable result,flag = 1
	//if(x >= off_x)
	//	result = NF_gauss(sig,off_x,x)
	//else
	//	result = -1.0 * NF_gauss(sig,off_x,x)
	//endif
	return -( 1+erf( (x-off_x)/ (sqrt(2*sig^2)) ) )/2 + 1
End
static Function NF_gauss(sig,off_x,x)
	variable sig,off_x,x
	return exp(-((x - off_x)^2)/(2*sig^2))/(sqrt(2*pi)*sig)
End

Function Noise_fil()
	//string /G NF_default_datafolder = "root:"
	DFREF tmp = GetNoiseFilter_Garbage()
	DFREF save_dfr = getDataFolderDFR()
	setdatafolder tmp
	string /G NF_Start_wave_name=""
	string /G  NF_End_wave_name=""
	variable /G NF_cut_width = 10
	variable /G NF_gauss_sigma = 100
	variable /G NF_manual_flag = 0
	setdatafolder save_dfr
	Noise_fil__()
End

static Function Noise_fil__()
	svar /SDFR=GetNoiseFilter_Garbage() NF_Start_wave_name,NF_End_wave_name
	nvar /SDFR=GetNoiseFilter_Garbage() NF_cut_width,NF_gauss_sigma
	DFREF now_save,dfr
	now_save = getdatafolderDFR()
	if (wintype("Noise_filter") != 0)
		DoWindow/F Noise_filter
	else
		NewPanel/W = (50,50,240,280)/K=1 as "Noise_filter"
		RenameWindow $S_name, Noise_filter
		
		GroupBox index_Group, pos={10,5},size={170,70},title="select waves"
		SetVariable setStartwave pos={20,25},size={150,20},value=NF_Start_wave_name,title="s_wave :"
		SetVariable setEndwave pos={20,49},size={150,20},value=NF_End_wave_name,title="e_wave :"
		string tmp_str = get_all_data_folder()
		//string trace_names = "\"" + wavelist("*",";","") + "\""
		//PopupMenu NF_popTrace_S, mode=1, Value=#trace_names,title="S_wave",pos={20,25},size={150,20}
		//PopupMenu NF_popTrace_E, mode=1, Value=#trace_names,title="E_wave",pos={20,49},size={150,20}
		string step_functions = "\"" + "step;erf" + "\""
		PopupMenu NF_select_step, mode=1,Value=#step_functions,title="select function",pos={10,72},size={17,20}
		SetVariable setwave_name pos={10,102},size={150,20},value=NF_cut_width,format="%d",title="cut width :"
		SetVariable set_gauss_sigma pos={10,122},size={150,20},value=NF_gauss_sigma,format="%d",title="sigma :"
		checkbox NF_checkbox_one pos={10,142},title="manual on / off "
		string tmp_folder = "\"" + tmp_str + "\""
		PopupMenu NF_select_folder,mode=1,Value=#tmp_folder,title="select folder",pos={10,166},size={17,20},proc=PopupMenuAction
		
		Button draw_button1, pos={10,196},size={170,20},proc=Fourier_func1,title="Cut Noise!", fSize=10
		//Button draw_button2, pos={10,136},size={170,20},proc=Fourier_func2,title="Insert Zero & Filter(Step2)", fSize=10
	endif

end

Function PopupMenuAction (ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum	// which item is currently selected (1-based)
	String popStr		// contents of current popup item as string
	
	//killwindow  $"Noise_filter"
	//svar NF_default_datafolder
	//NF_default_datafolder = popStr
	//Noise_fil__()
End

static Function /s get_all_data_folder()
	DFREF ro = root:
	String tmp = DataFolderdir(1,ro)
	String get_str
	splitstring /E="^FOLDERS:(.*)" tmp, get_str
	get_str = replacestring(",",get_str,";")
	print(get_str)
	return "root;" + remove_ca(get_str)
End
static function /s remove_ca(str)
	string str
	variable i
	for(i=strlen(str)-1; i>=0 && IsNewLine(str[i]);i-=1)
		;
	Endfor
	return str[0,i]
End
static function IsNewLine(c)
	string c
	return char2num(c) == 13 || char2num(c) == 10
End


Function Fourier_func1(ctrlName) : ButtonControl
	String ctrlName
	Main_Noise_filter_make_graph()
	print ctrlName
end

Function Main_Noise_filter_make_graph()
	svar /SDFR =GetNoiseFilter_Garbage() NF_Start_wave_name,NF_End_wave_name
	string get_str,reg_exp,input_str
	DFREF save_dfr,dfr,tmp_folder,spe_dfr
	ControlInfo /W = $("Noise_filter") NF_select_folder
	if(V_Value == 0)
		doalert 0,"no datafolder"
		return 0
	endif
	string NFF_default_datafolder = s_value
	save_dfr  = getdatafolderDFR()
	NFF_default_datafolder = get_datafolder_from_str(NFF_default_datafolder)
	dfr = $(NFF_default_datafolder)
	tmp_folder = GetNoiseFilter_Garbage()
	wave /Z /SDFR=dfr s_wave = $(NF_Start_wave_name)
	wave /Z /SDFR=dfr e_wave = $(NF_End_wave_name)
	if(waveexists(s_wave) == 0 || waveexists(e_wave) == 0)
		doalert 0,"incrrect data folder or no such a name of the wave"
		return 0
	endif
	reg_exp = NF_Start_wave_name + "(.*)" + NF_End_wave_name
	setdatafolder dfr
	input_str = wavelist("*",";","")
	setdatafolder save_dfr
	splitstring /E=reg_exp input_str, get_str
	if(stringmatch(NF_Start_wave_name,NF_End_wave_name))
		get_str = NF_Start_wave_name
	else
		get_str = NF_Start_wave_name + get_str + NF_End_wave_name
	endif
	print get_str//
	nvar /SDFR=GetNoiseFilter_Garbage() NF_manual_flag
	controlinfo /W = $("Noise_filter") NF_checkbox_one
	variable flag_of_manual = V_Value
	if(flag_of_manual == 1)
			//make graph and select cursol
		if(NF_manual_flag == 0)
			string element_t = stringfromlist(0,get_str)
			string gomi_t = element_t + "_FT_tmp"
			wave /Z /SDFR=dfr tmp_w_t = $(element_t)
			FFT /OUT=1 /DEST=tmp_folder:$(gomi_t) tmp_w_t
			wave /Z /SDFR=tmp_folder tmp_w_t_1 = $(gomi_t)
			NF_subwindow_func(tmp_w_t_1)
			NF_manual_flag = 1
			return 0
		endif
			//startIndex = pcsr(A)
	else
		NF_manual_flag = 0
	endif
	variable i
	for(i=0;i<itemsinlist(get_str);i+=1)
		string element=stringfromlist(i,get_str)
		wave /Z /SDFR=dfr tmp_wave = $(element)
		string gomi = element + "_FT"
		wave /Z /SDFR=tmp_folder FT_wave = $gomi
		if(waveexists(FT_wave))
			Really_Kill_Waves(FT_wave,tmp_folder)
		endif
		FFT /OUT=1 /DEST=tmp_folder:$(gomi) tmp_wave
		//reverse fourier transform
		variable startIndex,endIndex
		wave /Z /SDFR=tmp_folder FT_wave = $gomi
		variable length = numpnts(FT_wave)
		//wave /Z tmp_log = log(FT_wave)
		//Duplicate/O FT_wave tmp_folder:$"_tmp" /wave=wave_log
		//string s =  "make /O /N=" + num2str(length) + " tmp_folder:_tmp"
		//execute s
		make /O /N=(length) tmp_folder:$"_tmp" /wave=wave_log
		variable j,index_max
		for(j=0;j<length;j+=1)
			wave_log[j] = log(real(FT_wave[j])*real(FT_wave[j]) + imag(FT_wave[j]) * imag(FT_wave[j]))
		endfor
		variable max_value
		for(j=str2num(num2istr(length/4));j<length;j+=1)
			if(max_value < wave_log[j])
				max_value = wave_log[j]
				index_max = j
			endif
		endfor
		nvar /SDFR=GetNoiseFilter_Garbage() NF_cut_width,NF_gauss_sigma
		startIndex = index_max - NF_cut_width
		if(NF_manual_flag == 1)
			try
				startIndex = str2num(num2istr(pcsr(A,"tmp_graph")));AbortOnRTE
			catch
				variable cferror = GetRTError(1)	
				print GetErrMessage(cferror)
				doalert 0,"There is no cursols on the graph"
				NF_manual_flag = 0
				return 1
			endtry
			//print WinName(0,1)
		endif
		print startIndex
		endIndex = length
		controlinfo /W = $("Noise_Filter") NF_select_step
		string val = s_value
		if(cmpstr("step",val) == 0)
			make_zero(tmp_folder,gomi,startIndex)
		elseif(cmpstr("erf",val) == 0)
			cutting_data_by_step_func(tmp_folder,gomi,startIndex,NF_gauss_sigma)
		else
			doalert 0,"no filtering"
		endif
		spe_dfr = GetNoiseFilterFolder()
		IFFT /DEST=spe_dfr:$(element + "_NF") tmp_folder:$gomi
	endfor
	
end

Function NF_subwindow_func(w)
	wave /Z w
	if (wintype("tmp_graph") != 0)
		DoWindow/F tmp_graph
	else
		DFREF dfr = GetNoiseFilter_Garbage()
		DFREF save_dfr = getdatafolderDFR()
		display /N=tmp_graph /W=(200,100,580,320) w
		//display /W=(200,100,580,320) w 
		ModifyGraph log(left) = 1
		//string graph_name = Winname(0,1)
		NewPanel /K=1 /N=NF_Controls /W=(100,0,0,100)/HOST=tmp_graph /EXT=1 as "Controls"
		ModifyPanel /W=tmp_graph#NF_Controls,noEdit=1
		Button Fitting_Button,pos={25,40},size={50,20},title="Do Cut",proc=NF_subwindow
		DoWindow/F tmp_graph
		showinfo /W = tmp_graph
		setdatafolder dfr
		string /G NF_tmp_wave_name = nameofwave(w)
		setdatafolder save_dfr
	endif
End

Function NF_subwindow(ctrlName) :ButtonControl
	string ctrlName
	DFREF dfr = GetNoiseFilter_Garbage()	
	Main_Noise_filter_make_graph()
	nvar /SDFR=dfr NF_manual_flag
	NF_manual_flag = 0
	svar /SDFR=dfr NF_tmp_wave_name
	wave /Z /SDFR=dfr tmp_w = $(NF_tmp_wave_name)
	Really_Kill_Waves(tmp_w,dfr)
	if(waveexists(tmp_w))
		killwaves tmp_w
	endif
	killwindow tmp_graph
End

Function Really_Kill_Waves(w,dfr)
	wave w
	DFREF dfr
	DFREF saved
	saved = getdatafolderDFR()
	setdatafolder dfr
	if (waveExists(w))
		string name = nameofwave(w)
		removeWaveFromAllGraphs(name)
		removeWaveFromAllTables(name)
		killwaves w
	endif
	setdatafolder saved
End

Function removeWaveFromAllGraphs(name)
	string name
	string graphs = WinList("*",";","WIN:1")
	variable isDisplayed
	variable i
	for(i=0;i<itemsinlist(graphs);i+=1)
		string graph=stringfromlist(i,graphs)
		do
			checkDisplayed /w=$graph $name
				isDisplayed = V_flag
			if(isDisplayed)
				removeFromGraph /w=$graph $name
			endif
		while(isDisplayed)
	endfor
End

Function removeWaveFromAllTables(name)
	string name
	string tables=WinList("*",";","WIN:2")
	variable i
	for(i=0;i<itemsinlist(tables);i+=1)
		string table=stringfromlist(i,tables)
		checkDisplayed /w=$table $name
		if (V_flag == 1)
			RemoveFromTable /Z/W=$table $name
		endif
	endfor
End


static Function /DF GetNoiseFilterFolder()
	DFREF dfr = root:NoiseFilter
	if(DataFolderRefStatus(dfr) != 1)
		NewDataFolder /O root:NoiseFilter
		dfr = root:NoiseFilter
	endif
	return dfr
End

static Function /DF GetNoiseFilter_Garbage()
	DFREF dfr = root:NoiseFilter:Garbage
	if(Datafolderrefstatus(dfr)!=1)
		NewDataFolder /O root:NoiseFilter
		NewDataFolder /O root:NoiseFilter:Garbage
		dfr = root:NoiseFilter:Garbage
	endif
	return dfr
End


//--------------------------------------fitting and make absorbance----------------------------
//---new update plan---
//make button for merging spactra
//i need interface for merging spectra
Function /DF GetDFREF()
	DFREF dfr = root:Conv_Abs:For_Fitting
	if(DataFolderRefStatus(dfr) != 1)
		DFREF dfr = CreateFittingFolder()
	endif
	return dfr
End

static Function /DF CreateFittingFolder()
	NewDataFolder /O root:Conv_Abs
	NewDataFolder /O root:Conv_Abs:For_Fitting
	DFREF dfr = root:Conv_Abs:For_Fitting
	return dfr
End

static Function /DF Absorbance_Folder()
	DFREF dfr = root:Conv_Abs:spectra
	if(DataFolderRefStatus(dfr) != 1)
		dfr = root:Conv_Abs
		if(Datafolderrefstatus(dfr) != 1)
			NewDataFolder /O root:Conv_Abs
		endif
		NewDataFolder /O root:Conv_Abs:spectra
		dfr = root:Conv_Abs:spectra
	endif
	return dfr
End

static Function return_points(w,var)		//if var=0 then number of 1 bit, if var=1 then number of min_x, if var=2 then number of max_x
	wave /W w
	variable var
	if(waveexists(w)==0)
		print "wave does not exist"
		return -1
	endif
	variable length = numpnts(w)
	variable i, result
	switch(var)
		case 0:
			result = 0
			for(i=0;i<length;i+=1)
				if(w[i] == 1)
					result += 1
				endif
			endfor
			return result
		case 1:
			result = 0
			for(i=0;i<length;i+=1)
				if(w[i] == 1)
					result = i
					return result
				endif
			endfor
			return -1
		case 2:
			result=-1
			for(i=0;i<length;i+=1)
				if(w[i]==1)
					result = i
				endif
			endfor
			if(result == -1)
				return -1
			endif
			return result
		default:
			print "error of switch"
			return -1
	endswitch
End

Function Make_Absorbance()
	DFREF tmp =  GetDFREF()
	DFREF save_dfr = getdatafolderdfr()
	setdatafolder tmp
	variable /G global_execute_counter = 0
	variable /G show_counter=0
	string /G MA_save_func = "line"
	string /G MA_save_suffix = "0"
	setdatafolder save_dfr
	Make_Absorbance_()
End

static Function /S get_NF_save_func(func_names,str2)
	string func_names,str2
	variable i
	string result = str2
	for(i=0;i<itemsinlist(func_names);i+=1)
		string element = stringfromlist(i,func_names)
		if(cmpstr(element,str2)==0)
			continue
		endif
		result += ";" + element
	endfor
	return result
End
	

static Function Make_Absorbance_()
	if(strlen(Winlist("*",";","WIN:1"))==0)
		doalert 0,"Fitting requires a trace plotted in a graph window"
		return 0
	endif
	string graphStr = WinName(0,1)
	if (wintype(graphStr + "#Fitting_panel"))
		dowindow /F $(graphStr)	
	else
	//maybe this function needs clear procedure
	
	//make panel
	NewPanel /K=1 /N=Fitting_panel /W=(200,0,0,400)/HOST=$(graphStr)/EXT=1 as "Fitting Controls"
	ModifyPanel /W=$graphStr#Fitting_panel,noEdit=1
	variable i=0,deltaY=25,vL=30,vT=5,font=12,groupw=180
	
	GroupBox group0,pos={vL-20,vT+deltaY*i},size={groupw,deltaY*2},title="Data wave"
	GroupBox group0,fSize=font
	i+=1
	string tmp_str = BL_traces()
	tmp_str = "\"" + tmp_str + "\""
	PopupMenu popTrace, mode=1, Value=#tmp_str,title="",pos={vL,vT+deltaY*i},size={130,20}
	//PopupMenu popTrace, help={"select data wave"},proc=BL_popup
	i+=1.5
	GroupBox group1,pos={vL-20,vT+deltaY*i},size={groupw,deltaY*2},title="Baseline Type"
	GroupBox group1,fSize=font
	i+=1
	svar /SDFR=GetDFREF() MA_save_func, MA_save_suffix
	string fNames = "line;poly 3;poly 4;gauss;gauss3;lor;exp;dblexp;sin;hillequation;sigmoid;power;lognormal;"
	fNames = "\"" + get_NF_save_func(fNames,MA_save_func) + "\""
	PopupMenu popBL, mode=1, title="",pos={vL,vT+deltaY*i},size={130,20},Value=#fNames
	i+=1.5
	GroupBox group2,pos={vL-20,vT+deltaY*i},size={groupw,deltaY*2.5},title="Do Fit Operation"
	GroupBox group2,fSize=font
	i+=1
	Button Add_points,pos={vL-5,vT+deltaY*i+5},size={60,20},title="Add points",proc=HO_Add_points,fSize=10	//w_mask 0 -> 1
	Button Clear_all,pos={vL+70,vT+deltaY*i+5},size={70,20},title="Clear points",proc=HO_Clear_points,fSize=10 //w_mask = 0
	i+=2
	string save_name = "0;1;2;3;4;5"
	string save_name_t = "\"" + get_NF_save_func(save_name, MA_save_suffix) + "\""
	Popupmenu popsaveName,mode=1,Value=#save_name_t,title="",pos={vL,vT+deltaY*i},size={130,20}
	i+=1
	Button Fitting_Button,pos={60,vT+deltaY*i},size={70,20},title="Do Fit",proc=Fitting_Func
	i+=2
	DFREF tmp_ref = getdatafolderDFR()
	DFREF tmp = Absorbance_Folder()
	setdatafolder tmp
	// i wanna remove Freq ~ from the list
	string global_trace1 = wavelist("*",";","")
	variable j,k
	string result = ""
	string loop_str = "0;1;2;3;4;5"
	ControlInfo /W = $(graphStr + "#Fitting_panel") popTrace
	tmp_str = s_value
	for(j=0;j<itemsinlist(global_trace1);j+=1)
		string element=stringfromlist(j,global_trace1)
		string tmp1,tmp2
		splitstring /E="^(abs_of_)(.*)" element, tmp1,tmp2
		if(cmpstr("abs_of_",tmp1) == 0)
			for(k=0;k<itemsinlist(loop_str);k+=1)
				string element2 = stringfromlist(k,loop_str)
				if(cmpstr(tmp2,tmp_str + "_" + element2) == 0)
					result += element + ";"
				endif
			endfor
		endif
	endfor
	global_trace1 = "\"" + result + "\""
	setdatafolder tmp_ref
	PopupMenu popTrace1,mode=1, Value=#global_trace1,title="",pos={vL/4,vT+deltaY*i},size={130,20}
	i+=1
	PopupMenu popTrace2,mode=1,Value=#global_trace1,title="",pos={vL/4,vT+deltaY*i},size={130,20}
	i+=1
	Button merge,pos={vL-5,vT+deltaY*i+5},size={60,20},title="merge",proc=HO_Merge
	endif
	nvar /SDFR=GetDFREF() global_execute_counter
	HO_init(global_execute_counter)
End

Function HO_Merge(ctrlName):ButtonControl
	string ctrlName
	string str1,str2
	DFREF dfr = Absorbance_Folder()
	string graphStr=BL_getGraph()
	string panelStr = graphStr + "#Fitting_panel"
	controlinfo /W = $(panelStr) popTrace1
	if(V_Flag == 0)
		print "error"
		return 0
	endif
	str1 = s_value
	controlinfo /W = $(panelStr) popTrace2
	if(V_Flag == 0)
		print "error"
		return 0
	endif
	str2 = s_value
	controlinfo /W=$(panelStr) popTrace
	if(V_Flag == 0)
		doalert 0,"error in HO_merge function"
		return 0
	endif
	print s_value
	HO_Merge_(str1,str2,s_value,dfr)
	print ctrlName
End

Function HO_Add_points(ctrlName): ButtonControl
	String ctrlName
	DFREF dfr = GetDFREF()
	wave /Z /SDFR=dfr w_mask = w_mask
	if(waveexists(w_mask) == 0)
		return 0
	endif
	HO_make_w_mask(w_mask)
	print ctrlName
End

Function HO_Clear_points(ctrlName): ButtonControl
	string ctrlName
	DFREF dfr = GetDFREF()
	wave /Z /SDFR=dfr w_mask = w_mask
	if(waveexists(w_mask) == 0)
		return 0
	endif
	w_mask = 0
	nvar /SDFR=dfr show_counter
	wave /Z /SDFR=dfr w_base = w_base
	if(waveexists(w_base) && show_counter==1)
		show_counter = 0
		removefromgraph w_base
	endif
	print ctrlName
End

Function HO_init(num)
	variable num
	if(num ==1)
		return 0
	endif
	DFREF dfr = GetDFREF()
	string graphStr=BL_getGraph()
	string panelStr = graphStr + "#Fitting_panel"
	controlinfo /W = $(panelStr) popTrace
	if(V_Flag == 0)
		print "error you should append data in this graph"
		return 0
	Endif
	wave /Z w_data =TraceNameToWaveRef(graphStr,s_value)
	if (waveexists(w_data) == 0)
		print "none data"
		return 0
	endif
	
	duplicate /O w_data dfr:w_mask
	wave /Z /SDFR=dfr w_mask=w_mask
	w_mask = 0
End

static Function HO_make_w_mask(w)	// w[i]  value sets to 1  this function needs error handling
	Wave /Z w
	variable i
	if(waveexists(w) == 0)
		print "error in HO_make_w_mask_func"
		return 1
	endif
	try 
		i = pcsr(A);AbortOnRTE
	catch
		variable cferror = GetRTError(1)
		doalert 0,"there is no cursor in this graph"
		return 1
	endtry
	w[i] = 1
	return 0
End

Function Fitting_Func(ctrlName): ButtonControl
	String ctrlName
	string graphStr = BL_getGraph()
	string panelStr = graphStr+"#Fitting_panel"
	controlinfo /W=$(panelStr) popBL
	if(V_Flag==0)
		print "error in Fitting_func code 1"
		return 0
	endif
	string type = s_value
	svar /SDFR=GetDFREF() MA_save_func
 	MA_save_func = type
	controlinfo /W=$(panelStr) popTrace
	if(V_Flag == 0)
		print "error in Fitting_func code 2"
		return 0
	endif
	string data_wave_name = s_value
	DFREF dfr = GetDFREF()
	wave /Z w_data=TraceNameToWaveRef(graphStr,data_wave_name)
	wave /Z w_x=XWaveRefFromTrace(graphStr,data_wave_name)
	wave /Z /SDFR=dfr w_mask = w_mask
	if(waveexists(w_data) ==0 && waveexists(w_x)==0 && waveexists(w_mask))
		print("error")
		return 0
	endif
	if(return_points(w_mask,0) < 2)
		doalert 0, "You should put more than 2 cursors on the graph"
		return 0
	endif
	duplicate /O w_data dfr:w_base /Wave=w_base
	//variable result = 1
	//result = make_w_mask(w_mask)
	//if(result == 1)
	//	doalert 0,"cursor is not on the graph"
	//	return 0
	//endif
	Fit_wapper(w_data,w_x,w_mask,w_base,type)
	nvar /SDFR=dfr show_counter
	if(waveexists(w_base) && show_counter==0)
		show_counter =1
		appendtograph /W=$(graphStr) w_base vs w_x
		ModifyGraph /W=$(graphStr) rgb(w_base)=(54693,30000,65535)
	endif
	
	//cutting waves
	wave /Z w_base2_y = cutting_wave(w_base,dfr,"w_base2_y",w_mask)
	wave /Z w_base2_x = cutting_wave(w_x,dfr,"w_base2_x",w_mask)
	wave /Z w_base2_signal = cutting_wave(w_data,dfr,"w_base2_signal",w_mask)
	//convert to absorbance
	if(waveexists(w_base2_y) && waveexists(w_base2_x) && waveexists(w_base2_signal))
		wave /Z absorbance = convert_to_absorbance(w_base2_y,w_base2_signal,w_base2_x,dfr)
		if(waveexists(absorbance))
			DFREF t_dfr = Absorbance_Folder()
			//movewave cannot override so i have to killwaves
			string name1 = nameofwave(absorbance)
			string name2 = nameofwave(w_base2_x)
			wave /Z /SDFR=t_dfr tmp_wave1 = $(name1)
			wave /Z /SDFR=t_dfr tmp_wave2 = $(name2)
			DFREF saveDFR = getdatafolderDFR()
			setdatafolder t_dfr
			if(waveexists(tmp_wave1))
				Really_Kill_Waves(tmp_wave1,t_dfr)
			endif
			if(waveexists(tmp_wave2))
				Really_Kill_Waves(tmp_wave2,t_dfr)
			endif 
			
			movewave absorbance, t_dfr
			movewave w_base2_x, t_dfr
			controlinfo /W=$(panelStr) popsaveName
			string tmp_name1 =  "abs_of_" + data_wave_name + "_" + s_value
			string tmp_name2 = "Freq_of_" + data_wave_name + "_" + s_value
			svar /SDFR=GetDFREF() MA_save_suffix
			MA_save_suffix = s_value
			wave /Z tmp__1 = $(tmp_name1)
			wave /Z tmp__2 = $(tmp_name2)
			if(waveexists(tmp__1))
				Really_Kill_Waves(tmp__1,t_dfr)
			endif
			if(waveexists(tmp__2))
				Really_Kill_Waves(tmp__2,t_dfr)
			endif
			string tmp 
			tmp = "rename " + "w_base2_absorbance, " + tmp_name1
			execute tmp
			 tmp = "rename " + "w_base2_x, " + tmp_name2
			execute tmp
			setdatafolder saveDFR
			
		else
			print("Error making absorbance")
		endif 
	else
		print("ERROR")
	endif
	//display w_data vs w_x
	//print type
	//print data_wave_name
	setdatafolder saveDFR
 	KillWindow $(graphStr + "#Fitting_panel")
 	nvar /SDFR=dfr global_execute_counter
 	global_execute_counter = 1
 	Make_Absorbance_()
	print ctrlName
		
End
static Function /wave convert_to_absorbance(w_base,w_signal,w_Freq,dfr)
	wave /Z w_base,w_signal,w_Freq
	DFREF dfr
	variable i
	if(waveexists(w_base) && waveexists(w_signal) && waveexists(w_Freq))
		duplicate /O w_base dfr:w_base2_absorbance /wave=absorbance
		variable length = numpnts(w_base)
		for(i=0;i<length;i+=1)
			absorbance[i] = log(w_base[i]/w_signal[i])
		endfor
		return absorbance
	else
		print("error of converter")
	endif
	wave /Z t
	return t
End

static Function /wave cutting_wave(w,dfr,str,w_mask)
	wave /Z w
	DFREF dfr
	string str
	wave /Z w_mask
	variable index_s,index_e
	wave /Z tt
	if (waveexists(w))
		index_s = return_points(w_mask,1)
		index_e = return_points(w_mask,2)
		if(index_s != NaN && index_e != NaN)
			if(index_s > index_e)
				return tt
			else
				Duplicate /O /R=[index_s,index_e] w, dfr:$str /wave=tmp
				return tmp
			endif
		else
			return tt
		endif
	else
		print("error in cutting_wave")
		return tt
	endif
End
	
static Function make_w_mask(w)
	wave /Z w
	variable index_s,index_e
	try
		if (waveexists(w))
			index_s = pcsr(A);AbortOnRTE
			index_e = pcsr(B);AbortOnRTE
			if(index_s != NaN && index_e != NaN)
				if(index_s > index_e)
					doalert 0,"You should put A cursol forwarder than B cursol"
					return 1
				else
					w[index_s,index_e] = 1
					return 0
				endif
			else
				doalert 0,"You should put cursols on a graph!"
				print("error in make_w_mask")
				return 1
			endif
		else
			print("error in make_w_mask")
			return 1
		endif
	catch
		variable cferror = GetRTError(1)	
		print GetErrMessage(cferror)
		return 1
	endtry
	return 1
End
			
	
End	
	
static Function Fit_wapper(w,w_x,w_mask,w_out,fName)
	wave /Z w,w_x,w_mask,w_out
	string fName
	DFREF dfr = GetDFREF()
	DFREF save_dfr = getDatafolderDFR()
	setdatafolder dfr
	try
		strswitch(fName)
			case "line":
				CurveFit /Q/N line, w/M=w_mask /X=w_x /NWOK;AbortOnRTE
				wave w_coef
				w_out = waveexists(w_x) ? w_coef[0]+w_coef[1]*w_x : w_coef[0] + w_coef[1]*x
				break
			case "poly 3":
				CurveFit /Q/N poly 3, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				w_out = waveexists(w_x) ? poly(w_coef,w_x) : poly(w_coef,x)
				break
			case "poly 4":
				CurveFit /Q/N poly 4, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				w_out = waveexists(w_x) ? poly(w_coef,w_x) : poly(w_coef,x)
				break
			case "gauss":
				CurveFit /Q/N gauss, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				w_out = waveexists(w_x) ? gauss1D(w_coef,w_x) : gauss1D(w_coef,x)
				break
			case "gauss3":
				// first fit a gaussian to populate w_coef
				CurveFit /Q/N gauss, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				w_coef[0]=0 
				CurveFit /Q/N /H="1000" gauss, kwCWave=w_coef, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				w_out = waveexists(w_x) ? gauss1D(w_coef,w_x) : gauss1D(w_coef,x)
				break
			case "lor":
				CurveFit /Q/N lor, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				wave w=w_coef
				w_out = waveexists(w_x) ? w[0]+w[1]/((w_x-w[2])^2+w[3]) : w[0]+w[1]/((x-w[2])^2+w[3])
				break
			case "exp":
				CurveFit /Q/N exp, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				w_out = waveexists(w_x) ? w_coef[0]+w_coef[1]*exp(-w_coef[2]*w_x) : w_coef[0]+w_coef[1]*exp(-w_coef[2]*x)
				break
			case "dblexp":
				CurveFit /Q/N dblexp, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				wave w=w_coef
				w_out = waveexists(w_x) ? w[0]+w[1]*exp(-w[2]*w_x)+w[3]*exp(-w[4]*w_x) : w[0]+w[1]*exp(-w[2]*x)+w[3]*exp(-w[4]*x)
				break
			case "sin":
				CurveFit /Q/N sin, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				wave w=w_coef
				w_out = waveexists(w_x) ? w[0]+w[1]*sin(w[2]*w_x+w[3]) : w[0]+w[1]*sin(w[2]*x+w[3])
				break
			case "hillequation":
				CurveFit /Q/N hillequation, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				wave w=w_coef
				w_out = waveexists(w_x) ? w[0]+(w[1]-w[0])*(w_x^w[2]/(1+(w_x^w[2]+w[3]^w[2]))) : w[0]+(w[1]-w[0])*(x^w[2]/(1+(x^w[2]+w[3]^w[2])))
				break
			case "sigmoid":
				CurveFit /Q/N sigmoid, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				wave w=w_coef
				w_out = waveexists(w_x) ? w[0] + w[1]/(1+exp(-(w_x-w[2])/w[3])) : w[0] + w[1]/(1+exp(-(x-w[2])/w[3]))
				break
			case "power":
				CurveFit /Q/N power, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				w_out = waveexists(w_x) ? w_coef[0]+w_coef[1]*w_x^w_coef[2] : w_coef[0]+w_coef[1]*x^w_coef[2]
				break
			case "lognormal":
				CurveFit /Q/N lognormal, w /M=w_mask /X=w_x /NWOK; AbortOnRTE
				wave w_coef
				wave w=w_coef
				w_out = waveexists(w_x) ? w[0]+w[1]*exp(-(ln(w_x/w[2])/w[3])^2) : w[0]+w[1]*exp(-(ln(x/w[2])/w[3])^2)
				break
			default:
				break
		endswitch
	catch
		setdatafolder save_dfr
		if(V_AbortCode == -4)
			print "Error during curve fit:"
			variable cferror = getrterror(1)
			doalert 0,geterrmessage(cferror)
		endif
	endtry
	setdatafolder save_dfr
End

static Function /s BL_traces()
	string graphStr=BL_getGraph()
	string listStr=TraceNameList(graphStr,";",1+4)
	string removeStr="w_display;w_base;tangent0;tangent1;"
	removeStr += listmatch(listStr,"*_sub")
	removeStr += listmatch(listStr,"*_BL")
	listStr=RemoveFromList(removeStr,listStr,";",0)
	return listStr
End

static Function /T BL_getGraph()
	string listStr=WinList("*",";","WIN:1")
	string panelStr,graphStr
	variable i
	do
		graphStr = stringfromlist(i,listStr)
		panelStr = graphStr + "#Fitting_panel"
		if(wintype(panelStr))
			return graphStr
		endif
		i+=1
	while(i<itemsinlist(listStr)) 
	return ""
End

static Function HO_Merge_(str1,str2,outfile_name,dfr) //str1 ,str2 are wavenames
	string str1
	string str2
	string outfile_name
	DFREF dfr
	wave /Z /SDFR=dfr wave1_y = $(str1)
	wave /Z /SDFR=dfr wave2_y = $(str2)
	if(waveexists(wave1_y) == 0 || waveexists(wave2_y)==0)
		doalert 0, "input file y cannot open"
		return 0
	endif
	string str1_x_name,str2_x_name
	splitstring /E="^abs_of_(.*).*" str1, str1_x_name
	splitstring /E="^abs_of_(.*).*" str2, str2_x_name
	str1_x_name = "Freq_of_" + str1_x_name
	str2_x_name = "Freq_of_" + str2_x_name 
	wave /Z /SDFR=dfr wave1_x = $(str1_x_name)
	wave /Z /SDFR=dfr wave2_x = $(str2_x_name)
	if(waveexists(wave1_x) == 0 || waveexists(wave2_x)==0)
		doalert 0, "input file x cannot open"
		return 0
	endif
	
	variable wave1_x_number = numpnts(wave1_x)
	variable wave2_x_number = numpnts(wave2_x)
	//this function need to think about these details
	variable wave1_length = wave1_x[wave1_x_number-1] - wave1_x[0]
	variable wave2_length = wave2_x[wave2_x_number-1] - wave2_x[0]
	
	if(wave2_x[0] < wave1_x[0])
		HO_Merge_(str2,str1,outfile_name,dfr)
		return 0
	elseif(wave1_x[0] + wave1_length >= wave2_x[0])
		variable find_index = 0,i
		for(i=0;i<wave1_x_number;i+=1)
			if(wave2_x[0] == wave1_x[i])
				find_index = i
			endif
		endfor	
		variable total_len = find_index + wave2_x_number
		DFREF save_dfr = getDataFolderDFR()
		setdatafolder dfr
		make /O /N=(total_len) $("abs_of_" + outfile_name) /wave = output_w
		make /O /N=(total_len) $("Freq_of_" + outfile_name) /wave = out_x
		setdatafolder save_dfr
		for(i=0;i<find_index;i+=1)
			out_x[i] = wave1_x[i]
			output_w[i] = wave1_y[i]
		endfor
		for(i=find_index;i<total_len;i+=1)
			out_x[i] = wave2_x[i-find_index]
			output_w[i] = wave2_y[i-find_index]
		endfor
		if(cmpstr(str1,"abs_of_" + outfile_name + "_0") == 0)
			//Really_Kill_Waves(wave1_x,dfr)
			//Really_Kill_Waves(wave1_y,dfr)
		endif
		if(cmpstr(str2,"Freq_of_" + outfile_name + "_0") == 0)
			//Really_Kill_Waves(wave2_x,dfr)
			//Really_Kill_Waves(wave2_y,dfr)
		endif
		setdatafolder dfr
		//rename $("abs_of_" + outfile_name ) , $("abs_of_" + outfile_name + "_0")
		//rename $("$Freq_of_" + outfile_name) ,$("Freq_of_" + outfile_name + "_0")
		setdatafolder save_dfr
	else
		doalert 0,"these waves are not overwrapped"
	endif
	return 1
End

//-------------------Wrapper of NASU AutoFit Program----------------------
static Function /DF GetautoFitFolder()
	DFREF dfr = root:autoFit
	if(DataFolderRefStatus(dfr) != 1)
		NewDataFolder /O root:autoFit
		dfr = root:autoFit
	endif
	return dfr
End

static Function /DF Get_autoFitFolder_garbage()
	DFREF dfr = root:autoFit:garbage
	if(DataFolderRefstatus(dfr) != 1)
		dfr = root:autoFit
		if(DataFolderRefstatus(dfr) != 1)
		NewDatafolder /O root:autoFit
		endif
		NewDatafolder /O root:autoFit:garbage
		dfr = root:autoFit:garbage
	endif
	return dfr
End


static Function autoFit_Panel__()
	if (wintype("AutoFit_Panel_p") != 0)
		DoWindow/F AutoFit_Panel_p
	else
		NewPanel/W = (50,50,240,240)/K=1 as "AutoFit_Panel"
		RenameWindow $S_name, AutoFit_Panel_p
		DFREF save_dfr,dfr
		dfr = Get_autoFitFolder_garbage()
		save_dfr = getdatafolderDFR()
		setdatafolder dfr
		string /G AF_Start_wave_name=""
		string /G AF_End_wave_name=""
		variable /G AF_num_of_peaks = 0
		setdatafolder save_dfr
		GroupBox index_Group, pos={10,5},size={170,70},title="select waves"
		SetVariable setStartwave pos={20,25},size={150,20},value=AF_Start_wave_name,title="s_wave :"
		SetVariable setEndwave pos={20,49},size={150,20},value=AF_End_wave_name,title="e_wave :"
		string tmp_str = get_all_data_folder()
		checkbox NS_autoFit_checkbox pos={10,82},title="fix on / off (need const freq)"
		SetVariable setwave_name pos={10,102},size={150,20},value=AF_num_of_peaks,format="%d",title="num of peaks :"
		string tmp_folder = "\"" + tmp_str + "\""
		PopupMenu NF_select_folder,mode=1,Value=#tmp_folder,title="select waves folder",pos={10,126},size={17,20},proc=PopupMenuAction_for_autofit
		Button draw_button_NS_AutoFit, pos={10,156},size={170,20},proc=NS_Auto_Fit_button,title="Do Fit", fSize=10
	endif
End
static Function Main_autoFit_Panel()
	DFREF out_dfr = GetautoFitFolder()
	DFREF garbage_dfr = Get_autoFitFolder_garbage()
	svar /SDFR=garbage_dfr AF_Start_wave_name,AF_End_wave_name
	controlinfo /W = $("AutoFit_Panel_p") NF_select_folder
	string spec_folder = get_datafolder_from_str(s_value)
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox
	variable flag_for_fix = V_Value
	variable index_s,index_e
	DFREF dfr = $(spec_folder)
	DFREF save_dfr = getDataFolderDFR()
	nvar /SDFR=garbage_dfr AF_num_of_peaks
	AF_make_W_coef(garbage_dfr,save_dfr,AF_num_of_peaks)
	
	
	wave /Z /SDFR=garbage_dfr w_coef = $("W_coef")
	if(waveexists(w_coef) == 0)
		doalert 0,"you should try to choose the folder which is W_coef"
		return 1
	endif
	
	svar /SDFR=garbage_dfr AF_Start_wave_name,AF_End_wave_name
	wave /Z /SDFR=dfr s_wave = $(AF_Start_wave_name)
	wave /Z /SDFR=dfr e_wave = $(AF_End_wave_name)
	if(waveexists(s_wave) == 0 || waveexists(e_wave) == 0)
		doalert 0,"incrrect data folder or no such a name of the wave"
		return 0
	endif
	string get_str
	string reg_exp = AF_Start_wave_name + "(.*)" + AF_End_wave_name
	setdatafolder dfr
	string input_str = wavelist("*",";","")
	setdatafolder save_dfr
	splitstring /E=reg_exp input_str, get_str
	if(stringmatch(AF_Start_wave_name,AF_End_wave_name))
		get_str = AF_Start_wave_name
	else
		get_str = AF_Start_wave_name + get_str + AF_End_wave_name
	endif
	
	if(AF_num_of_peaks < 1 || AF_num_of_peaks > 10)
		doalert 0,"error of peaks; You should try to reduce num of peaks or increase ones"
		return 1
	endif
	
	try
		index_s = pcsr(A);AbortOnRTE
		index_e = pcsr(B);AbortOnRTE
		if(index_s != NaN && index_e != NaN)
			if(index_s > index_e)
				doalert 0,"You should put A cursol forwarder than B cursol"
				return 1
			else
				if(flag_for_fix == 1)
					autoFit(index_s,index_e,w_coef,dfr,out_dfr,AF_num_of_peaks,get_str,fix=1)
				else
					autoFit(index_s,index_e,w_coef,dfr,out_dfr,AF_num_of_peaks,get_str)
				endif
			endif
		else
			doalert 0,"You should put cursols on a graph!"
			print("error in make_w_mask")
			return 1
		endif
		
	catch
		variable cferror = GetRTError(1)
		doalert 0,"You should put cursols on a graph!"
		print GetErrMessage(cferror)
		return 1
	endtry
End
Function autoFit_Panel()
	autoFit_Panel__()
End

Function NS_Auto_Fit_button(ctrlName) : ButtonControl
	string ctrlName
	Main_autoFit_Panel()
End

Function PopupMenuAction_for_autofit(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum	// which item is currently selected (1-based)
	String popStr		// contents of current popup item as string
	//killwindow  $"AutoFit_Panel_p"
	//svar NF_default_datafolder
	//NF_default_datafolder = popStr
	//Noise_fil__()
	//autoFit_Panel()
End

static Function /S dfr2str(dfr)
	DFREF dfr
	String dfstr
	DFREF initialFolder = GetDataFolderDFR()
	SetDataFolder dfr
	dfstr = GetDataFolder(1)
	SetDataFolder initialFolder
	return dfstr
End

static Function MyFit(n,x0,x1,W_coef,drawFlag,fix,fix2,length,dfr,out_dfr,freq_name,spec_name)
	variable n,x0,x1,drawFlag,fix,length
	string fix2,freq_name,spec_name
	DFREF dfr,out_dfr
	wave /Z W_coef
	string s, H
	if(waveexists(W_coef) == 0)
		doalert 0,"cannot read wave W_coef"
		return 0
	endif
	DFREF garbage_dfr = Get_autoFitFolder_garbage()
	wave /Z /SDFR=dfr Freq = $(freq_name)
	//w_sigma is current folder
	DFREF save_dfr = getdatafolderDFR()
	wave /Z /SDFR=save_dfr W_sigma = $("W_sigma")
	
	duplicate /O/D dfr:$(spec_name) ,garbage_dfr:$("wt")
	wave /Z /SDFR=garbage_dfr w = $("wt")
	if(waveexists(w) == 0 || waveexists(Freq) == 0 || waveexists(W_sigma) == 0)
		doalert 0,"cannot read wt or freq or w_sigma"
		return 0
	endif
	WaveStats /Q W_coef
	variable size = DimSize(W_coef,0)
	string tmp_str,tmp_str1="\""
	setdatafolder garbage_dfr
	if (fix==1)
		tmp_str = "FuncFit/H=" + tmp_str1 + NS_GET_BIT(size) + tmp_str1+"/Q/NTHR=0 "+NS_GET_STR(size)+" " +dfr2str(garbage_dfr)+nameofwave(W_coef)+ " " +dfr2str(garbage_dfr)+nameofwave(w)+"[" + num2str(x0)+ ","+num2str(x1)+"]" +" /X=" + dfr2str(dfr) + "Freq /D" 
		execute tmp_str
	else
		tmp_str = "FuncFit/H="+ tmp_str1 + fix2 + tmp_str1 + "/Q/NTHR=0 " + NS_GET_STR(size) + " " + dfr2str(garbage_dfr) + nameofwave(W_coef) + " " +  dfr2str(garbage_dfr) + nameofwave(w) + "[" + num2str(x0) + "," + num2str(x1) + "]" +" /X=" + dfr2str(dfr) + "Freq /D" 
		execute tmp_str
	endif
	setdatafolder save_dfr
	wave /Z /SDFR=garbage_dfr fit_wt = $("fit_wt")
	if(waveexists(fit_wt) == 0)
		doalert 0,"not found fit_wt wave"
		return 0
	endif
	string tmp1 = "Integral_",tmp2 = "Int_sigma_",tmp3 = "Freq_"
	if  ( size >= 4 )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp1 + num2str(0))
		tmp_wave[n] = W_coef[1]*W_coef[3]*sqrt(pi)  //+ Int2
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp2 + num2str(0))
		tmp_wave[n] = sqrt( W_coef[1]*W_coef[1]*pi*W_sigma[3]*W_sigma[3] + W_coef[3]*W_coef[3]*pi*W_sigma[1]*W_sigma[1] )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp3 + num2str(0))
		tmp_wave[n] = W_coef[2] 
	endif
	if ( size >= 8 )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp1 + num2str(1))
		tmp_wave[n] = W_coef[5]*W_coef[7]*sqrt(pi)  //+ Int2
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp2 + num2str(1))
		tmp_wave[n] = sqrt( W_coef[5]*W_coef[5]*pi*W_sigma[7]*W_sigma[7] + W_coef[7]*W_coef[7]*pi*W_sigma[5]*W_sigma[5] )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp3 + num2str(1))
		tmp_wave[n] = W_coef[6] 
	endif
	if ( size >= 11 )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp1 + num2str(2))
		tmp_wave[n] = W_coef[8]*W_coef[10]*sqrt(pi)  //+ Int2
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp2 + num2str(2))
		tmp_wave[n] = sqrt( W_coef[8]*W_coef[8]*pi*W_sigma[10]*W_sigma[10] + W_coef[10]*W_coef[10]*pi*W_sigma[8]*W_sigma[8] )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp3 + num2str(2))
		tmp_wave[n] = W_coef[9] 
	endif
	if ( size >= 14 )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp1 + num2str(3))
		tmp_wave[n] = W_coef[11]*W_coef[13]*sqrt(pi)  //+ Int2
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp2 + num2str(3))
		tmp_wave[n] = sqrt( W_coef[11]*W_coef[11]*pi*W_sigma[13]*W_sigma[13] + W_coef[13]*W_coef[13]*pi*W_sigma[11]*W_sigma[11] )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp3 + num2str(3))
		tmp_wave[n] = W_coef[12] 
	endif
	if ( size >= 17 )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp1 + num2str(4))
		tmp_wave[n] = W_coef[14]*W_coef[16]*sqrt(pi)  //+ Int2
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp2 + num2str(4))
		tmp_wave[n] = sqrt( W_coef[14]*W_coef[14]*pi*W_sigma[16]*W_sigma[16] + W_coef[16]*W_coef[16]*pi*W_sigma[14]*W_sigma[14] )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp3 + num2str(4))
		tmp_wave[n] = W_coef[15] 
	endif
	if ( size >= 20 )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp1 + num2str(5))
		tmp_wave[n] = W_coef[17]*W_coef[19]*sqrt(pi)  //+ Int2
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp2 + num2str(5))
		tmp_wave[n] = sqrt( W_coef[17]*W_coef[17]*pi*W_sigma[19]*W_sigma[19] + W_coef[19]*W_coef[19]*pi*W_sigma[17]*W_sigma[17] )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp3 + num2str(5))
		tmp_wave[n] = W_coef[18] 
	endif
	if ( size >= 23 )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp1 + num2str(6))
		tmp_wave[n] = W_coef[20]*W_coef[22]*sqrt(pi)  //+ Int2
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp2 + num2str(6))
		tmp_wave[n] = sqrt( W_coef[20]*W_coef[20]*pi*W_sigma[22]*W_sigma[22] + W_coef[22]*W_coef[22]*pi*W_sigma[20]*W_sigma[20] )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp3 + num2str(6))
		tmp_wave[n] = W_coef[21] 
	endif
	if ( size >= 26 )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp1 + num2str(7))
		tmp_wave[n] = W_coef[23]*W_coef[25]*sqrt(pi)  //+ Int2
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp2 + num2str(7))
		tmp_wave[n] = sqrt( W_coef[23]*W_coef[23]*pi*W_sigma[25]*W_sigma[25] + W_coef[25]*W_coef[25]*pi*W_sigma[23]*W_sigma[23] )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp3 + num2str(7))
		tmp_wave[n] = W_coef[24] 
	endif
	if ( size >= 29 )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp1 + num2str(8))
		tmp_wave[n] = W_coef[26]*W_coef[28]*sqrt(pi)  //+ Int2
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp2 + num2str(8))
		tmp_wave[n] = sqrt( W_coef[26]*W_coef[26]*pi*W_sigma[28]*W_sigma[28] + W_coef[28]*W_coef[28]*pi*W_sigma[26]*W_sigma[26] )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp3 + num2str(8))
		tmp_wave[n] = W_coef[27] 
	endif
	if ( size >= 32 )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp1 + num2str(9))
		tmp_wave[n] = W_coef[29]*W_coef[31]*sqrt(pi)  //+ Int2
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp2 + num2str(9))
		tmp_wave[n] = sqrt( W_coef[29]*W_coef[29]*pi*W_sigma[31]*W_sigma[31] + W_coef[31]*W_coef[31]*pi*W_sigma[29]*W_sigma[29] )
		wave /Z /SDFR=out_dfr tmp_wave = $(tmp3 + num2str(9))
		tmp_wave[n] = W_coef[30] 
	endif
	wave /Z /SDFR=out_dfr tmp_wave = $(tmp1 + num2str(0))
	print "[" + spec_name + "]" + " y0="  + num2str(W_coef[0]) + ", A=" + num2str(W_coef[1]) + ", x0=" + num2str(W_coef[2]) + ", width=" + num2str(W_coef[3]) + ", k=" + num2str(W_coef[4]) + ", S1=" + num2str(tmp_wave[n]) 
	
	if(drawFlag == 1)
		Display fit_wt; AppendToGraph w[x0,x1] vs Freq[x0,x1]
	endif
End

	
static Function /S NS_GET_BIT(i)
	variable i
	switch(i)
		case 4:
			return "0010"
		case 5:
			return "00100"
		case 8:
			return "00100010"
		case 11:
			return "00100010010"
		case 14:
			return "00100010010010"
		case 17:
			return "00100010010010010"
		case 20:
			return "00100010010010010010"
		case 23:
			return "00100010010010010010010"
		case 26:
			return "00100010010010010010010010"
		case 29:
			return "00100010010010010010010010010"
		case 32:
			return "00100010010010010010010010010010"
		default:
			doalert 0,"error of Fitting You should reduce waves"
			return ""
	endswitch
End
static Function /s NS_GET_STR(i)
	variable i
	switch(i)
		case 4:
			return "gauss1"
		case 5:
			return "gauss_linear"
		case 8:
			return "gauss2"
		case 11:
			return "gauss3"
		case 14:
			return "gauss4"
		case 17:
			return "gauss5"
		case 20:
			return "gauss6"
		case 23:
			return "gauss7"
		case 26:
			return "gauss8"
		case 29:
			return "gauss9"
		case 32:
			return "gauss10"
		default:
			doalert 0,"error of Fitting You should reduce waves"
			return ""
	endswitch
End

Function autoFit(x0,x1,W_coef,dfr,out_dfr,length,wave_names,[rev,fix,subtract,fix2]) //length is the number of peaks
	variable x0,x1,rev,fix,length		//dfr is where wt is  
	DFREF dfr,out_dfr
	wave W_coef,subtract
	string fix2,wave_names
	variable i,j,isSubtract = 1
	if ( ParamIsDefault(fix) )
		fix = 0
	endif
	if ( ParamIsDefault(rev) )
		rev = 0
	endif
	if ( ParamIsDefault(subtract) )
		isSubtract = 0
	endif
	if ( ParamIsDefault(fix2) )
		fix2 = ""
	endif
	variable len_of_waves = 0
	for(i=0;i<itemsinlist(wave_names);i+=1)
		len_of_waves += 1
	endfor
	for(i=0;i<length;i+=1)
		string tmp,tmp1= "Integral_;Int_sigma_;Freq_"
		for(j=0;j<itemsinlist(tmp1);j+=1)
			string element=stringfromlist(j,tmp1)
			tmp = element + num2str(i)
			make /O /D /N=(len_of_waves) out_dfr:$(tmp)
			//SetScale x i0,i1+1, out_dir:$(tmp)
		endfor
	endfor
	string merge_name = ""
	if(rev == 0)
		;
	else
		//wave_names sort reversely
		for(i=0;i<itemsinlist(wave_names);i+=1)
			string element1 = stringfromlist(i,wave_names)
			merge_name = element1 +";"+ merge_name
		endfor
		wave_names = merge_name
	endif
	for(i=0;i<itemsinlist(wave_names);i+=1)
		string element2=stringfromlist(i,wave_names)
		MyFit(i,x0,x1,W_coef,0,fix,fix2,length,dfr,out_dfr,"Freq",element2)
		DFREF garbage_dfr = Get_autoFitFolder_garbage()
		wave /Z /SDFR=garbage_dfr editWave = $("wt")
		if(waveexists(editWave) == 0)
			doalert 0,"wt wave cannot open"
			return 0
		endif
		if(isSubtract == 1)
			for(j=0;j<dimsize(subtract,0);j+=1)
				editWave[j] -= (W_coef[1]*exp(-((subtract(j) - W_coef[2])/W_coef[3])^2))
			endfor
			duplicate /O editWave,dfr:$(element2)
		endif
	endfor
End

Function gauss_linear(w,x):FitFunc
	Wave w
	variable x
	return w[0]+w[1]*exp(-( (x-w[2])/w[3] )^2 ) + w[4]*x
End
Function gauss1(w,x) : FitFunc
	Wave w
	Variable x
	return w[0]+w[1]*exp(-( (x-w[2])/w[3] )^2 )
End
Function gauss2(w,x) : FitFunc
	Wave w
	Variable x
	return w[0]+w[1]*exp(-( (x-w[2])/w[3] )^2 ) +w[5]*exp(-( (x-w[6])/w[7] )^2 )  + w[4]*x
End
Function gauss3(w,x) : FitFunc
	Wave w
	Variable x
	return w[0]+w[1]*exp(-( (x-w[2])/w[3] )^2 ) +w[5]*exp(-( (x-w[6])/w[7] )^2 )  +w[8]*exp(-( (x-w[9])/w[10] )^2 ) + w[4]*x
End
Function gauss4(w,x) : FitFunc
	Wave w
	Variable x
	return w[0]+w[1]*exp(-( (x-w[2])/w[3] )^2 ) +w[5]*exp(-( (x-w[6])/w[7] )^2 )  +w[8]*exp(-( (x-w[9])/w[10] )^2 )  +w[11]*exp(-( (x-w[12])/w[13] )^2 ) + w[4]*x
End
Function gauss5(w,x) : FitFunc
	Wave w
	Variable x
	return w[0]+w[1]*exp(-( (x-w[2])/w[3] )^2 ) +w[5]*exp(-( (x-w[6])/w[7] )^2 )  +w[8]*exp(-( (x-w[9])/w[10] )^2 )  +w[11]*exp(-( (x-w[12])/w[13] )^2 ) + w[14]*exp(-( (x-w[15])/w[16] )^2 ) +w[4]*x
End
Function gauss6(w,x) : FitFunc
	Wave w
	Variable x
	return w[0]+w[1]*exp(-( (x-w[2])/w[3] )^2 ) +w[5]*exp(-( (x-w[6])/w[7] )^2 )  +w[8]*exp(-( (x-w[9])/w[10] )^2 )  +w[11]*exp(-( (x-w[12])/w[13] )^2 ) + w[14]*exp(-( (x-w[15])/w[16] )^2 )+ w[17]*exp(-( (x-w[18])/w[19] )^2 ) +w[4]*x
End
Function gauss7(w,x) : FitFunc
	Wave w
	Variable x
	return w[0]+w[1]*exp(-( (x-w[2])/w[3] )^2 ) +w[5]*exp(-( (x-w[6])/w[7] )^2 )  +w[8]*exp(-( (x-w[9])/w[10] )^2 )  +w[11]*exp(-( (x-w[12])/w[13] )^2 ) + w[14]*exp(-( (x-w[15])/w[16] )^2 )+ w[17]*exp(-( (x-w[18])/w[19] )^2 ) + w[20]*exp(-( (x-w[21])/w[22] )^2 ) +w[4]*x
End
Function gauss8(w,x) : FitFunc
	Wave w
	Variable x
	return w[0]+w[1]*exp(-( (x-w[2])/w[3] )^2 ) +w[5]*exp(-( (x-w[6])/w[7] )^2 )  +w[8]*exp(-( (x-w[9])/w[10] )^2 )  +w[11]*exp(-( (x-w[12])/w[13] )^2 ) + w[14]*exp(-( (x-w[15])/w[16] )^2 )+ w[17]*exp(-( (x-w[18])/w[19] )^2 ) + w[20]*exp(-( (x-w[21])/w[22] )^2 ) +w[23]*exp(-( (x-w[24])/w[25] )^2 ) +w[4]*x
End
Function gauss9(w,x) : FitFunc
	Wave w
	Variable x
	return w[0]+w[1]*exp(-( (x-w[2])/w[3] )^2 ) +w[5]*exp(-( (x-w[6])/w[7] )^2 )  +w[8]*exp(-( (x-w[9])/w[10] )^2 )  +w[11]*exp(-( (x-w[12])/w[13] )^2 ) + w[14]*exp(-( (x-w[15])/w[16] )^2 )+ w[17]*exp(-( (x-w[18])/w[19] )^2 ) + w[20]*exp(-( (x-w[21])/w[22] )^2 ) +w[23]*exp(-( (x-w[24])/w[25] )^2 ) +w[26]*exp(-( (x-w[27])/w[28] )^2 ) +w[4]*x
End
Function gauss10(w,x) : FitFunc
	Wave w
	Variable x
	return w[0]+w[1]*exp(-( (x-w[2])/w[3] )^2 ) +w[5]*exp(-( (x-w[6])/w[7] )^2 )  +w[8]*exp(-( (x-w[9])/w[10] )^2 )  +w[11]*exp(-( (x-w[12])/w[13] )^2 ) + w[14]*exp(-( (x-w[15])/w[16] )^2 )+ w[17]*exp(-( (x-w[18])/w[19] )^2 ) + w[20]*exp(-( (x-w[21])/w[22] )^2 ) +w[23]*exp(-( (x-w[24])/w[25] )^2 ) +w[26]*exp(-( (x-w[27])/w[28] )^2)+w[29]*exp(-( (x-w[30])/w[31] )^2)+w[4]*x
End

static Function AF_make_W_coef(dfr,peak_fit_dfr,size)
	DFREF dfr,peak_fit_dfr
	variable size
	DFREF save_dfr = getdatafolderDFR()
	variable tmp1 = 2 + size * 3
	wave /Z /SDFR=peak_fit_dfr base_w = $("Baseline Coefs")
	if(waveexists(base_w) == 0)
		doalert 0,"baseline does not exist"
		return 0
	endif
	if(dimsize(base_w,0) != 2)
		doalert 0,"You need to select Linear in Baseline"
		return 0
	endif
	variable i
	wave /Z /T /SDFR=peak_fit_dfr tmp_w = $("MPF2_ResultsListWave")
	if(waveexists(tmp_w) == 0)
		doalert 0,"Please press Peak Results button"
		return 0
	endif
	string result = ""
	for(i=0;i<size;i+=1)
		if(i==0)
			result += num2str(base_w[0])
		endif
		//wave /Z /T t_tmp_w = tmp_w[i]
		string flag = tmp_w[i][1]
		if(cmpstr(flag,"Gauss") != 0)
		 	doalert 0,"Please select Gauss in Peak Type"
		 	return 0
		endif
		result += "," + tmp_w[i][14] 
		result += "," + tmp_w[i][10]
		result += "," + tmp_w[i][12]
		if(i==0)
			result += "," + num2str(base_w[1])
		endif
		 
	endfor
	
	setdatafolder dfr
	string tmp = "make /O/D/N=(" + num2istr(tmp1) + ") W_coef = {" + result + "}"
	execute tmp
	//print tmp
	setdatafolder save_dfr
End

//---------------subtract------------

static Function /DF Subtract_Folder()
	DFREF dfr = root:Subtract
	if(DataFolderRefStatus(dfr) != 1)
		NewDataFolder /O root:Subtract
		dfr = root:Subtract
	endif
	return dfr
End

Function subtract()
	if(strlen(Winlist("*",";","WIN:1"))==0)
		doalert 0,"Subtraction requires at least two traces in the graph window"
		return 0
	endif
	string graphStr = WinName(0,1)
	if(wintype(graphStr + "#Subtract_panel"))
		dowindow /F $(graphStr)
	else
		NewPanel /K=1 /N=Subtract_panel /W=(200,0,0,200) /HOST=$(graphStr) /EXT=1 as "Subtract Controls"
		ModifyPanel /W=$graphStr#Subtract_panel,noEdit=1
		variable i=0,deltaY=25,vL=30,vT=5,font=12,groupw=180
		GroupBox group3,pos={vL-20,vT+deltaY*i},size={groupw,deltaY*2},title="Select waves"
		GroupBox group3,fSize=font
		i+=1
		string tmp_str = "\"" + SB_traces() + "\""
		PopupMenu SB_popTrace1,mode=1,Value=#tmp_str,title="",pos={vL,vT+deltaY*i},size={130,20}
		i+=2
		PopupMenu SB_popTrace2,mode=1,Value=#tmp_str,title="",pos={vL,vT+deltaY*i},size={130,20}
		i+=2
		Button Fitting_Button,pos={60,vT+deltaY*i},size={90,20},title="Do Subtract",proc=SB_main
	endif
End

static Function /s SB_traces()
	string graphStr=WinName(0,1)
	string listStr=TraceNameList(graphStr,";",1+4)
	string removeStr="w_display;w_base;tangent0;tangent1;"
	removeStr += listmatch(listStr,"*_sub")
	removeStr += listmatch(listStr,"*_BL")
	listStr=RemoveFromList(removeStr,listStr,";",0)
	return listStr
End

Function SB_main(ctrlName)
	string ctrlName
	string graphStr = WinName(0,1)
	ControlInfo /W = $(graphStr + "#Subtract_panel") SB_popTrace1
	string wave1_name = s_value
	ControlInfo /W = $(graphStr + "#Subtract_panel") SB_popTrace2
	string wave2_name = s_value
	if(cmpstr(wave1_name,wave2_name)==0)
		doalert 0, "same traces"
		return 0
	endif
	wave /Z wave1_y=TraceNameToWaveRef(graphStr,wave1_name)
	wave /Z wave1_x=XWaveRefFromTrace(graphStr,wave1_name)
	wave /Z wave2_y=TraceNameToWaveRef(graphStr,wave2_name)
	wave /Z wave2_x=XWaveRefFromTrace(graphStr,wave2_name)
	variable i,tmp_index=-1
	variable flag = 0
	DFREF saved_dfr,out_dfr = Subtract_Folder()
	saved_dfr = getdatafolderDFR()
	if(wave1_x[0] < wave2_x[0])
		for(i=0;i<numpnts(wave1_x);i+=1)
			if(wave1_x[i] == wave2_x[0])
				tmp_index = i
				flag = 1
				break
			endif
		endfor
	elseif(wave1_x[0] == wave2_x[0])
		tmp_index = 0
	else
		for(i=0;i<numpnts(wave2_x);i+=1)
			if(wave2_x[i] == wave1_x[0])
				tmp_index = i
				flag = -1
				break
			endif
		endfor
	endif
	if(tmp_index==-1)
		doalert 0,"cannot find overwrapping"
		return 0
	endif
	variable max_index
	switch(flag)
		case 0:
					variable len1 = numpnts(wave1_x)
					variable len2 = numpnts(wave2_x)
					variable length 
					if(len1 > len2)
						length = len2
					else
						length = len1
					endif
					string new_wave = wave1_name +"_sub"
					setdatafolder out_dfr
					make /O /N=(length) $(new_wave + "_freq") /wave=out_wave_x
					make /O /N=(length) $(new_wave + "_abs") /wave=out_wave_y
					setdatafolder saved_dfr
					if(len1 > len2)
						for(i=0;i<len2;i+=1)
							out_wave_y[i] = wave1_y[i] - wave2_y[i]
							out_wave_x[i] = wave2_x[i]
						endfor
					else
						for(i=0;i<len1;i+=1)
							out_wave_y[i] = wave1_y[i] - wave2_y[i]
							out_wave_x[i] = wave1_x[i]
						endfor
					endif
					break
			case 1:
					variable len3 = numpnts(wave1_x)
					variable len4 = numpnts(wave2_x)
					variable result_index = -1
					for(i=0;i<len3;i+=1)
						if(wave1_x[i] == wave2_x[0])
							result_index = i
						endif
					endfor
					if(result_index == -1)
						doalert 0,"not overwrapped"
						return 0
					endif
					variable len3_tmp = len3 - result_index
					string new_wave2 = wave1_name +"_sub"
					variable length2
					if(len3_tmp > len4)
						length2 = len4
					else
						length2 = len3_tmp
					endif
					setdatafolder out_dfr
					make /O /N=(length2) $(new_wave2 + "_freq") /wave=out_wave_x
					make /O /N=(length2) $(new_wave2+ "_abs") /wave=out_wave_y
					setdatafolder saved_dfr
					if(len3_tmp > len4)
						for(i=result_index;i<len4+result_index;i+=1)
							out_wave_x[i-result_index] = wave2_x[i - result_index]
							out_wave_y[i-result_index] = wave1_y[i] - wave2_y[i-result_index]
						endfor
					else
						for(i=0;i<len3_tmp;i+=1)
							out_wave_x[i] = wave2_x[i]
							out_wave_y[i] = wave1_y[i+result_index] - wave2_y[i]
						endfor
					endif
					break
			case -1:
					print "--------"
					variable len5 = numpnts(wave1_x)
					variable len6 = numpnts(wave2_x)
					variable result_index2 = -1
					for(i=0;i<len6;i+=1)
						if(wave2_x[i] == wave1_x[0])
							result_index2 = i
						endif
					endfor
					if(result_index2 == -1)
						doalert 0,"not overwrapped  error code = -1"
						return 0
					endif
					variable len6_tmp = len6 - result_index2
					string new_wave3 = wave1_name +"_sub" 
					variable length3 
					if(len6_tmp > len5)
						 length3 = len5
					else
						length3 = len6_tmp
					endif
					print len5 - result_index2
					print length3
					setdatafolder out_dfr
					make /O /N=(length3) $(new_wave3 + "_freq") /wave=out_wave_x
					make /O /N=(length3) $(new_wave3+ "_abs") /wave=out_wave_y
					setdatafolder saved_dfr
					if(len6_tmp > len5)
						for(i=result_index2;i<len5+result_index2;i+=1)
							out_wave_x[i-result_index2] = wave1_x[i-result_index2]
							out_wave_y[i-result_index2] = wave1_y[i-result_index2] - wave2_y[i]
						endfor
					else
						for(i=0;i<len6_tmp;i+=1)
							out_wave_x[i] = wave1_x[i]
							out_wave_y[i] = wave1_y[i] - wave2_y[i+result_index2]
						endfor
					endif					
	endswitch
	appendtograph /W=$graphStr out_wave_y vs out_wave_x
	print wave1_name
	print wave2_name
End
	

Function test()
	//print ChildWindowList("graph0")
	//Dowindow /F $("graph0#name0")
	//string tmp = wavelist("*",";","WIN:")
	DFREF ro = root:
	string tmp = datafolderdir(1,ro)
	print tmp
End

Function test1()
	wave /Z w_coef = root:$("W_coef")
	autoFit(pcsr(A),pcsr(B),w_coef,root:,root:Untitled,3,"jws180402004;jws180402005;jws180402006",fix=1)
end

Function test3()
	print pcsr(A,"tmp_graph")
End

