#pragma rtGlobals=3	// Use modern global access method and strict wave access.
#include <WaveSelectorWidget>
#include <PopupWaveSelector>
#include <Multi-peak Fitting 2.0>
#include <HierarchicalListWidget>

#if (IgorVersion() < 7) 
	Menu "Noise_Filter"
	"Noise_fil"
	"Noise_fil_2"
	"Make_Absorbance"
	"subtract"
	"CUT_DATA"
	"Create_calculated_wave"
	"Graph_Offset_Easily"
	"Correlation"
	"autoFit_Panel"
	"New_AutoFit_Function"
	End
#else
	Menu "Noise_Filter"
	"Noise_fil"
	"Noise_fil_2"
	"Make_Absorbance"
	"subtract"
	"CUT_DATA"
	"Create_calculated_wave"
	"Graph_Offset_Easily"
	"Correlation"
	"autoFit_Panel_(Only_supported_on_Igor6)"
	"New_AutoFit_Function_(Only_supported_on_Igor6)"
	End
	
#endif


//-------------------------renewal NF(changed select waves)-----------------

Function Noise_fil_2()
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
	Noise_fil_2_core()
End

Function NF_PopupWaveSelectorNotify(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	DFREF now_save = getDataFolderDFR()
	setdatafolder GetNoiseFilter_Garbage()
	svar NF_Start_wave_name,NF_End_wave_name
	setdatafolder now_save
	NF_Start_wave_name = wavepath
	NF_End_wave_name = wavepath
	print NF_Start_wave_name
	//print "Selected wave:",wavepath, " using control", ctrlName
end

static Function Noise_fil_2_core()
	DFREF now_save,dfr
	now_save = getdatafolderDFR()
	setdatafolder GetNoiseFilter_Garbage()
	svar NF_Start_wave_name,NF_End_wave_name
	nvar NF_cut_width,NF_gauss_sigma
	setdatafolder now_save
	if(wintype("Noise_filter") != 0)
		killwindow Noise_filter
	endif
	if (wintype("Noise_filter2") != 0)
		//DoWindow/F Noise_filter2
		
		Killwindow Noise_filter2
		Noise_fil_2()
		return -2
	else
		NewPanel/W = (50,50,240,280)/K=1 as "Noise_filter2"
		RenameWindow $S_name, Noise_filter2
		
		GroupBox index_Group, pos={10,5},size={170,70},title="select waves"
		//SetVariable setStartwave pos={20,25},size={150,20},value=NF_Start_wave_name,title="s_wave :"
		//SetVariable PopupWaveSelectorSV3,pos={20,25},size={150,70},title="s_wave"
		//SetVariable PopupWaveSelectorSV3,bodyWidth= 150
		//MakeSetVarIntoWSPopupButton("Noise_filter2", "PopupWaveSelectorSV3", "NF_PopupWaveSelectorNotify", "root:Packages:DemoPopupWaveSelector:WSSetVar3")
		Button NF2_Button,pos={20,25},size={150,20}
		MakeButtonIntoWSPopupButton("Noise_filter2", "NF2_Button", "NF_PopupWaveSelectorNotify", options=PopupWS_OptionFloat)
		string tmp_str = get_all_data_folder()
		string step_functions = "\"" + "step;erf" + "\""
		PopupMenu NF_select_step, mode=1,Value=#step_functions,title="select function",pos={10,72},size={17,20}
		SetVariable setwave_name pos={10,102},size={150,20},value=NF_cut_width,format="%d",title="cut width :"
		SetVariable set_gauss_sigma pos={10,122},size={150,20},value=NF_gauss_sigma,format="%d",title="sigma :"
		checkbox NF_checkbox_one pos={10,142},title="manual on / off "
		string tmp_folder = "\"" + tmp_str + "\""
		PopupMenu NF_select_folder,mode=1,Value=#tmp_folder,title="select folder",pos={10,166},size={17,20},proc=PopupMenuAction
		
		Button draw_button1, pos={10,196},size={170,20},proc=Fourier_func2,title="Cut Noise!", fSize=10
		
	endif

end


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
	if(wintype("Noise_filter2") != 0)
		killwindow Noise_filter2
	endif
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

Function Fourier_func2(ctrlName):ButtonControl
	String ctrlName
	Main_Noise_filter_make_graph2("Noise_filter2")
	print ctrlName
End

Function Fourier_func1(ctrlName) : ButtonControl
	String ctrlName
	Main_Noise_filter_make_graph("Noise_filter")
	print ctrlName
end

Function Main_Noise_filter_make_graph2(my_panel_name)
	string my_panel_name
	svar /SDFR =GetNoiseFilter_Garbage() NF_Start_wave_name,NF_End_wave_name
	string get_str,reg_exp,input_str
	DFREF save_dfr,tmp_folder,spe_dfr
	save_dfr  = getdatafolderDFR()
	tmp_folder = GetNoiseFilter_Garbage()
	wave   s_wave = $(NF_Start_wave_name)
	wave   e_wave = $(NF_End_wave_name)
	if(waveexists(s_wave) == 0 || waveexists(e_wave) == 0)
		doalert 0,"incrrect data folder or no such a name of the wave"
		return 0
	endif
	if(mod(numpnts(s_wave) ,2) != 0)
		doalert 0, "length of wave you selected must be even."
		return 0
	endif
	reg_exp = NF_Start_wave_name + "(.*)" + NF_End_wave_name
	//splitstring /E=reg_exp input_str, get_str
	//print get_str//
	nvar /SDFR=GetNoiseFilter_Garbage() NF_manual_flag
	controlinfo /W = $(my_panel_name) NF_checkbox_one
	variable flag_of_manual = V_Value
	if(flag_of_manual == 1)
			//make graph and select cursol
		if(NF_manual_flag == 0)
			string element_t = NF_Start_wave_name
			wave ele_w = $(NF_Start_wave_name)
			string gomi_t = nameofwave(ele_w) + "_FT_tmp"
			wave tmp_w_t = ele_w
			
			FFT /OUT=1 /DEST=tmp_folder:$(gomi_t) tmp_w_t 
			wave /Z /SDFR=tmp_folder tmp_w_t_1 = $(gomi_t)
			
			NF_subwindow_func2(tmp_w_t_1)
			NF_manual_flag = 1
			return 0
		endif
			//startIndex = pcsr(A)
	else
		NF_manual_flag = 0
	endif

		string element=NF_Start_wave_name
		wave  tmp_wave = $(NF_Start_wave_name)
		string gomi = nameofwave(tmp_wave)+ "_FT"
		wave /Z /SDFR=tmp_folder FT_wave = $gomi
		if(waveexists(FT_wave))
			Really_Kill_Waves(FT_wave,tmp_folder)
		endif
		FFT /OUT=1 /DEST=tmp_folder:$(gomi) tmp_wave
		//reverse fourier transform
		variable startIndex,endIndex
		wave /Z /SDFR=tmp_folder FT_wave = $gomi
		variable length = numpnts(FT_wave)
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
				startIndex = str2num(num2istr(pcsr(A,"tmp_graph2")));AbortOnRTE
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
		controlinfo /W = $(my_panel_name) NF_select_step
		string val = s_value
		if(cmpstr("step",val) == 0)
			make_zero(tmp_folder,gomi,startIndex)
		elseif(cmpstr("erf",val) == 0)
			cutting_data_by_step_func(tmp_folder,gomi,startIndex,NF_gauss_sigma)
		else
			doalert 0,"no filtering"
		endif
		spe_dfr = GetNoiseFilterFolder()
		print nameofwave(tmp_wave)
		IFFT /DEST=spe_dfr:$(nameofwave(tmp_wave) + "_NF") tmp_folder:$gomi
	
	
end


Function Main_Noise_filter_make_graph(my_panel_name)
	string my_panel_name
	svar /SDFR =GetNoiseFilter_Garbage() NF_Start_wave_name,NF_End_wave_name
	string get_str,reg_exp,input_str
	DFREF save_dfr,dfr,tmp_folder,spe_dfr
	ControlInfo /W = $(my_panel_name) NF_select_folder
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
	controlinfo /W = $(my_panel_name) NF_checkbox_one
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
		controlinfo /W = $(my_panel_name) NF_select_step
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

Function NF_subwindow_func2(w)
	wave /Z w
	if (wintype("tmp_graph2") != 0)
		DoWindow/F tmp_graph2
	else
		DFREF dfr = GetNoiseFilter_Garbage()
		DFREF save_dfr = getdatafolderDFR()
		string showIng_wave = NameOfWave(w) + "_s"
		setDataFolder dfr
		make /O /N=(numpnts(w)) $(showIng_wave) /wave=show_wa
		setDataFolder save_dfr
		variable i
		for(i=0;i<numpnts(w);i+=1)
			show_wa[i] = sqrt(real(w[i])^2 + imag(w[i])^2)
		endfor
		display /N=tmp_graph2 /W=(200,100,580,320) show_wa
		//display /W=(200,100,580,320) w 
		ModifyGraph log(left) = 1
		//string graph_name = Winname(0,1)
		NewPanel /K=1 /N=NF_Controls /W=(100,0,0,100)/HOST=tmp_graph2 /EXT=1 as "Controls"
		ModifyPanel /W=tmp_graph2#NF_Controls,noEdit=1
		Button Fitting_Button,pos={25,40},size={50,20},title="Do Cut",proc=NF_subwindow2
		DoWindow/F tmp_graph2
		showinfo /W = tmp_graph2
		setdatafolder dfr
		string /G NF_tmp_wave_name = nameofwave(w)
		setdatafolder save_dfr
	endif
End

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

Function NF_subwindow2(ctrlName) :ButtonControl
	string ctrlName
	DFREF dfr = GetNoiseFilter_Garbage()	
	Main_Noise_filter_make_graph2("Noise_filter2")	//-------------------------------------直さなければならない
	nvar /SDFR=dfr NF_manual_flag
	NF_manual_flag = 0
	svar /SDFR=dfr NF_tmp_wave_name
	wave /SDFR=dfr tmp_w = $(NF_tmp_wave_name)
	Really_Kill_Waves(tmp_w,dfr)
	if(waveexists(tmp_w))
		killwaves tmp_w
	endif
	killwindow tmp_graph2
End

Function NF_subwindow(ctrlName) :ButtonControl
	string ctrlName
	DFREF dfr = GetNoiseFilter_Garbage()	
	Main_Noise_filter_make_graph("Noise_filter")	//-------------------------------------直さなければならない
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
	


// ----------------------------cut data----------------------------------------------------

Function CUT_DATA()	//公開用関数
	HO_cut_data()
End

static Function HO_cut_data()
	string graphStr = WinName(0,1)
	if(wintype(graphStr + "#Cutting_data_panel"))
		dowindow /F $(graphStr)
	else
		NewPanel /K=1 /N=Cutting_data_panel /W=(180,0,0,155) /HOST=$(graphStr) /EXT=1 as "Cutting Controls"
		ModifyPanel /W=$graphStr#Cutting_data_panel,noEdit=1
		variable i=0,deltaY=25,vL=30,vT=5,font=12,groupw=180
		GroupBox group3,pos={vL-20,vT+deltaY*i},size={groupw,deltaY*2},title="put A and B cursors "
		GroupBox group3,fSize=font
		i+=1
		Button Fitting_Button,pos={45,vT+deltaY*i},size={90,20},title="Do Cut",proc=HO_cut_data_core
		i+=1
		checkbox HO_cut_checkbox pos={10,vT+deltaY*i},title="yデータ毎にxデータ保存?"
		i+=1
		GroupBox group4,pos={vL-30,vT+deltaY*i},size={groupw,deltaY*2},title="説明 : このマクロはグラフに描かれている"
		i+=1
		GroupBox group5,pos={vL-30,vT+deltaY*i},size={groupw,deltaY*2},title="全てのデータに対し A,B間でカットし、"
		i+=1
		GroupBox group6,pos={vL-30,vT+deltaY*i},size={groupw,deltaY*2},title="Cut_dataフォルダにそのデータが入る"
	endif
End

Function HO_cut_data_core(ctrlName)
	string ctrlName
	string graphStr = WinName(0,1)
	variable a_i,b_i,flag = 0
	try
		a_i = pcsr(A,graphStr);AbortOnRTE
	catch
		variable cferror = GetRTError(1)
		flag = 1
	endtry
	try
		b_i = pcsr(B,graphStr);AbortOnRTE
	catch
		cferror = GetRTError(1)
		if(flag == 0)
			flag = 2
		else
			flag = 3
		endif
	endtry
	switch(flag)
		case 1:
			doalert 0,"There isn't A cursor in this graph(" + graphStr + ")"
			return 1
		case 2:
			doalert 0,"There isn't B cursor in this graph(" + graphStr + ")"
			return 1
		case 3:
			doalert 0,"There aren't A & B cursors in this graph(" + graphStr + ")"
			return 1
		default:
			break
	endswitch
	variable num_of_points = b_i - a_i
	if(num_of_points < 0)
		variable tmp_num
		tmp_num = a_i
		a_i = b_i
		b_i = tmp_num
		num_of_points = b_i - a_i
	elseif(num_of_points == 0)
		doalert 0,"AとBの間隔がありません"
		return 1
	endif
	controlinfo /W=$(graphStr + "#Cutting_data_panel") HO_cut_checkbox
	variable flag_of_each = V_Value
	string listStr = TraceNameList(graphStr,";",1+4)
	variable i
	for(i=0;i<itemsinlist(listStr);i+=1)		//ただしX軸がある場合全部同じデータ数でなければならない
		string wave_name = stringfromlist(i,listStr)
		wave /Z tmp_y = TraceNameToWaveRef(graphStr,wave_name)
		wave /Z tmp_x = XwaveRefFromTrace(graphStr,wave_name)
		//print(num2str(num_of_points) + "," + num2str(numpnts(tmp_y)))
		if(numpnts(tmp_y) < num_of_points+2)
			doalert 0, "This wave("  + nameofwave(tmp_y) + ")doesn't have enough data to cut!!\n ,But the others are Okay"	//カットするデータの方が小さい時はエラーを出す
			continue
		endif
		if(waveexists(tmp_x) == 0)	//x軸データはない
			//その実装
			DFREF default_folder = GetDataFolderDFR()
			string x_wave_name_1 = wave_name + "_x"
			setDataFolder cut_data_Folder()
			if(flag_of_each != 0 || i == 0)
				make /O /N=(num_of_points + 1) $(x_wave_name_1)	//x
			endif
			make /O /N=(num_of_points + 1) $(wave_name)			// y
			setDataFolder default_folder
			if(flag_of_each != 0  || i == 0)
				wave /SDFR=cut_data_Folder() out_x = $(x_wave_name_1)
			endif
			wave /SDFR=cut_data_Folder() out_y = $(wave_name)
			variable l,m
			for(l=a_i,m=0;l<=b_i;l+=1,m+=1)
				if(flag_of_each != 0 || i == 0)
					out_x[m] = l
				endif
				out_y[m] = tmp_y[l]
			endfor
		else
			DFREF default_folder = GetDataFolderDFR()
			//string x_wave_name = nameofwave(tmp_x)
			string x_wave_name = wave_name + "_x"
			SetDataFolder cut_data_Folder()
			make /O /N = (num_of_points + 1) $wave_name
			if(flag_of_each != 0 || i == 0)
				make /O /N = (num_of_points + 1) $(x_wave_name)
			endif
			SetDataFolder default_folder
			if(flag_of_each != 0  || i == 0)
				wave /SDFR=cut_data_Folder() output_x = $(x_wave_name)
			endif
			wave /SDFR=cut_data_Folder() output_y = $wave_name
			variable j,k
			for(j=a_i,k=0;j<=b_i;j+=1,k+=1)
				if(flag_of_each != 0 || i == 0)
					output_x[k] = tmp_x[j]
				endif
				output_y[k] = tmp_y[j]
			endfor
		endif
	endfor
End

static Function /DF cut_data_Folder()
	DFREF dfr = root:Cut_data
	if(DataFolderRefStatus(dfr) != 1)
		NewDataFolder /O root:Cut_data
		dfr = root:Cut_data
	endif
	return dfr
End

//========================= resize & capture window size & save =====================
//アスペクト比をどうするか？

static Function get_window_size_width(graph_name)
	string graph_name
	getwindow /Z $(graph_name) psize
	variable result = V_right - V_left
	return result
End

static Function get_window_size_height(graph_name)
	string graph_name
	getwindow /Z $(graph_name) psize
	variable result = V_bottom - V_top
	return result
End

Function resize_window(graph_name,i_width,i_height)
	string graph_name
	variable i_width,i_height
	variable now_width = get_window_size_width(graph_name)
	variable calc_percentage = i_width / now_width * 100
	ModifyGraph /W=$(graph_name) gfMult=100
	ModifyGraph /W=$(graph_name) width=i_width
	ModifyGraph /W=$(graph_name) height=i_height
	
End



//-------------------------------------autofit written by hoitaku------------------------

// ----------------------------- Create _calculated_ wave ------------------------

static Function /DF get_calc_datafolder()		//Working Directory for Dem Functions or Demo variables
	DFREF dfr = root:calculated_wave_
	if(DataFolderRefStatus(dfr) != 1)
		NewDataFolder /O root:calculated_wave_
		dfr = root:calculated_wave_
	endif
	return dfr
End

Function CCW_PopupWaveSelector(event,wavepath,windowName,ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	DFREF now_save = getDataFolderDFR()
	setdatafolder get_calc_datafolder()
	svar selected_wave_name
	setdatafolder now_save
	selected_wave_name = wavepath
End

Function CCW_core(ctrlName):ButtonControl
	string ctrlName
	DFREF now_save = getDataFolderDFR()
	setdatafolder get_calc_datafolder()
	svar selected_wave_name
	setdatafolder now_save
	wave /Z source_wave = $(selected_wave_name)
	if(waveexists(source_wave)==0)
		doalert 0,"This wave ("  + nameofwave(source_wave) + ") doesn't exist"
		return -1
	endif
	print wavetype(source_wave,0)
	if((wavetype(source_wave) & 0x01)== 1)
		doalert 0,"This wave is complex wave, So this procedure can create the wave.\nbut don't guarantee correct axis wave"
	endif
	string output_name = nameofwave(source_wave) + "_calc_"		// wave name wasn't allowed taking more than 40 characters
	variable i,len
	len = numpnts(source_wave)
	setdatafolder get_calc_datafolder()
	make /O /N=(len) $(output_name)
	setdatafolder now_save
	wave /SDFR=get_calc_datafolder() out = $(output_name)
	variable start_offset = DimOffset(source_wave,0)
	variable delta_x = DimDelta(source_wave,0)
	for(i=0;i<len;i+=1)
		if(i==0)
			out[0] = start_offset
			continue
		endif
		out[i] = out[i-1] + delta_x
	endfor
	
End

Function Create_calculated_wave()
	DFREF current_dfr = GetDataFolderDFR()
	DFREF working_dfr = get_calc_datafolder()
	string Window_name = "Create_calculated"
	string button_name = "Execute_CCW"
	if(wintype(Window_name) != 0)
		Dowindow /F $(Window_name)
	else
		setdatafolder working_dfr
		string /G selected_wave_name = ""
		setdatafolder current_dfr
		NewPanel /W=(50,50,240,150) /K=1 as Window_name
		RenameWindow $(S_name), $(Window_name)
		GroupBox group1, pos={10,5},size={170,50},title="Select wave"
		Button $(button_name),pos={20,25},size={150,20}
		MakeButtonIntoWSPopupButton(Window_name,button_name,"CCW_PopupWaveSelector",options=PopupWS_OptionFloat)
		Button execute_button,pos={20,60},size={150,20},proc=CCW_core,title="Create",fSize=10
	endif
	
End

// ------------------------ Graph Offset Easily ------------------------

static Function /DF Get_Graph_Offset_DFR()
	DFREF dfr = root:Graph_Offset_Easily
	if(DataFolderRefStatus(dfr) != 1)
		NewDataFolder /O root:Graph_Offset_Easily
		dfr = root:Graph_Offset_Easily
	endif
	return dfr
End

static Function /s Graph_Offset_Traces(window_name)
	string window_name
	string listStr=TraceNameList(window_name,";",1+4)
	string removeStr = "w_display;w_base;tangent0;tangent1;"
	listStr = RemoveFromList(removeStr,listStr,";",0)
	return listStr
End

Function GO_Refresh_func(ctrlName):ButtonControl
	string ctrlName
	string graph_window_name = WinName(0,1)
	string SubPanel_name = "#Graph_Offset"
	svar /SDFR=Get_Graph_Offset_DFR() popup_str
	popup_str = "\"" + Graph_Offset_Traces(graph_window_name)+ "\""
	SetActiveSubwindow $(graph_window_name + SubPanel_name)
	PopupMenu GGO_popupTrace1 value=#(popup_str)
	ControlUpdate /W=$(graph_window_name + SubPanel_name) GGO_popupTrace1
End

Function Graph_Offset_Easily_core(ctrlName):ButtonControl
	string ctrlName
	string GGO_Window_Name = WinName(0,1)
	nvar /SDFR=Get_Graph_Offset_DFR() GGO_Start_Offset_X,GGO_Start_Offset_delta_X,GGO_Start_Offset_Y,GGO_Start_Offset_delta_Y,GGO_Last_Selected_wave_count
	svar /SDFR=Get_Graph_Offset_DFR() GGO_Last_Selected_wave
	variable x_offset = GGO_Start_Offset_X + GGO_Start_Offset_delta_X * GGO_Last_Selected_wave_count
	variable y_offset = GGO_Start_Offset_Y + GGO_Start_Offset_delta_Y * GGO_Last_Selected_wave_count
	string SubPanel_name = "#Graph_Offset"
	ControlInfo /W=$(GGO_Window_Name + SubPanel_name) GGO_popupTrace1
	string selected_wave = s_value
	//Print "--------" + selected_wave
	ModifyGraph /W=$(GGO_Window_Name) offset($(selected_wave))={x_offset,y_offset}
	controlinfo /W=$(GGO_Window_Name + SubPanel_name) set_auto_count_index
	variable flag_auto_count = V_Value
	//Print flag_auto_count
	if(flag_auto_count == 1)
		GGO_Last_Selected_wave_count += 1
	endif
	Doupdate
End

Function Graph_Offset_Easily()
	if(strlen(Winlist("*",";","WIN:1")) == 0)
		doalert 0,"This procedure needs at least one trace on this graph window"
		return -1
	endif
	string graph_window_name = WinName(0,1)
	string SubPanel_name = "#Graph_Offset"
	string SubPanel_name_main = "Graph_Offset"
	DFREF current_dfr = GetDatafolderDFR()
	
	nvar /SDFR=Get_Graph_Offset_DFR() /Z GGO_Start_Offset_X
	if(!Nvar_Exists(GGO_Start_Offset_X))
		setdatafolder Get_Graph_Offset_DFR()
		variable /G GGO_Start_Offset_X=0
		setdatafolder current_dfr
	endif
	nvar /SDFR=Get_Graph_Offset_DFR() /Z GGO_Start_Offset_delta_X
	if(!Nvar_Exists(GGO_Start_Offset_delta_X))
		setdatafolder Get_Graph_Offset_DFR()
		variable /G GGO_Start_Offset_delta_X=0
		setdatafolder current_dfr
	endif
	nvar /SDFR=Get_Graph_Offset_DFR() /Z GGO_Start_Offset_Y
	if(!Nvar_Exists(GGO_Start_Offset_Y))
		setdatafolder Get_Graph_Offset_DFR()
		variable /G GGO_Start_Offset_Y=0
		setdatafolder current_dfr
	endif
	nvar /SDFR=Get_Graph_Offset_DFR() /Z GGO_Start_Offset_delta_Y
	if(!Nvar_Exists(GGO_Start_Offset_delta_Y))
		setdatafolder Get_Graph_Offset_DFR()
		variable /G GGO_Start_Offset_delta_Y=0
		setdatafolder current_dfr
	endif
	svar /SDFR=Get_Graph_Offset_DFR() /Z GGO_Last_Selected_wave		//wave name	(in fact , this variable hasn't been used yet)
	if(!Svar_Exists(GGO_Last_Selected_wave))
		setdatafolder Get_Graph_Offset_DFR()
		string /G GGO_Last_Selected_wave="d"
		setdatafolder current_dfr
	endif
	nvar /SDFR=Get_Graph_Offset_DFR() /Z GGO_Last_Selected_wave_count		//前回と異なっていたらカウントしていく
	if(!Nvar_Exists(GGO_Last_Selected_wave_count))								//
		setdatafolder Get_Graph_Offset_DFR()
		variable /G GGO_Last_Selected_wave_count=0
		setdatafolder current_dfr
	endif
	
	svar /SDFR=Get_Graph_Offset_DFR() /Z popup_str
	if(!svar_exists(popup_str))
		setdatafolder Get_Graph_Offset_DFR()
		string /G popup_str	//For popup selector
		setdatafolder current_dfr
	endif
	
	if(wintype(graph_window_name + SubPanel_name))
		Dowindow /F $(graph_window_name)
		popup_str = "\"" + Graph_Offset_Traces(graph_window_name)+ "\""
		SetActiveSubwindow $(graph_window_name + SubPanel_name)
		//print ControlNameList(graph_window_name + SubPanel_name, ";", "*_popupTrace1")
		//Print "----- "
		//ModiFyControl GGO_popupTrace1 value=#popup_str
		//SetVariable GGO_popupTrace1 value=#popup_str
		PopupMenu GGO_popupTrace1 value=#(popup_str)
		ControlUpdate /W=$(graph_window_name + SubPanel_name) GGO_popupTrace1
		//expect Update later 
	else
		NewPanel /K=1 /N=$(SubPanel_name_main) /W=(200,0,0,280) /Host=$(graph_window_name) /EXT=1 as "Graph Offset"
		ModifyPanel /W=$(graph_window_name + SubPanel_name),noEdit=1
		variable i=0,deltaY=25,vL=30,vT=5,groupw=180
		Groupbox group1,pos={vL-20,vT+deltaY*i},size={groupw,deltaY*2},title="Select Wave"
		i+=1
		popup_str = "\"" + Graph_Offset_Traces(graph_window_name)+ "\""
		PopupMenu GGO_popupTrace1,mode=1,Value=#popup_str,title="",pos={vL,vT+deltaY*i},size={130,20}
		i+=1.2
		SetVariable setStart_Offset_X pos={vL-20,vT+deltaY*i},size={150,20},value=GGO_Start_Offset_X,title="Start X Offset:"
		i += 1
		SetVariable setStart_Offset_delta_X pos={vL-20,vT+deltaY*i},size={150,20},value=GGO_Start_Offset_delta_X,title="delta X            :"
		i += 1
		SetVariable setStart_Offset_Y pos={vL-20,vT+deltaY*i},size={150,20},value=GGO_Start_Offset_Y,title="Start Y Offset:"
		i += 1
		SetVariable setStart_Offset_delta_Y pos={vL-20,vT+deltaY*i},size={150,20},value=GGO_Start_Offset_delta_Y,title="delta Y            :"
		i+=1
		SetVariable set_index_Offset pos={vL-20,vT+deltaY*i},size={150,20},value=GGO_Last_Selected_wave_count,title="index              :"
		i+=1
		CheckBox set_auto_count_index pos={vL-20,vT+deltaY*i},title="auto increment count index?"
		i+=1
		Button Refresh_button,pos={60,vT+deltaY*i},size={90,20},title="Refresh Popup",proc=GO_Refresh_func
		i+=1
		Button Execute_Offset,pos={60,vT+deltaY*i},size={90,20},title="Offset !",proc=Graph_Offset_Easily_core
	endif
End	

// ------------------------------- Correlation -----------------------------
static Function /DF get_Correlation_datafolder()
	DFREF dfr = root:Correlation
	if(DataFolderRefStatus(dfr) != 1)
		NewDataFolder /O root:Correlation
		dfr = root:Correlation
	endif
	return dfr
End

Function CL_PopupWaveSelector_Bid(event,wavepath,windowName,ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	DFREF now_save = getDataFolderDFR()
	setdatafolder get_Correlation_datafolder()
	svar DT_Global_WaveName_Bid
	setdatafolder now_save
	DT_Global_WaveName_Bid = wavepath
End

Function CL_PopupWaveSelector_Ask(event,wavepath,windowName,ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	DFREF now_save = getDataFolderDFR()
	setdatafolder get_Correlation_datafolder()
	svar DT_Global_WaveName_Ask
	setdatafolder now_save
	DT_Global_WaveName_Ask = wavepath
End

Function Correlation_core(ctrlName) : ButtonControl
	string ctrlName
	DFREF now_dfr = getDataFolderDFR()
	DFREF working_dfr = get_Correlation_datafolder()
	string Window_name = "Demo_Trader"
	setdatafolder working_dfr
	svar DT_Global_WaveName_Bid,DT_Global_WaveName_Ask
	setdatafolder now_dfr
	wave /Z bid_wave = $(DT_Global_WaveName_Bid)
	wave /Z ask_wave = $(DT_Global_WaveName_Ask)
	if(waveexists(bid_wave) == 0 || waveexists(ask_wave) == 0)
		doalert 0,"cannot find " + "waves you selected" 
		return -1
	endif
	if(numpnts(bid_wave) != numpnts(ask_wave))
		doalert 0,"Different Length Data you selected\nThe Length of these waves must be same"
	endif
	variable i
	variable sum1_bid=0,sum2_bid=0,sum1_ask=0,sum2_ask=0
	for(i=0;i<numpnts(bid_wave);i+=1)
		sum1_bid += bid_wave[i]
		sum2_bid += bid_wave[i] * bid_wave[i]
		sum1_ask += ask_wave[i]
		sum2_ask += ask_wave[i] * ask_wave[i]
	endfor
	variable mean1_bid,mean2_bid,mean1_ask,mean2_ask
	mean1_bid = sum1_bid/numpnts(bid_wave)
	mean2_bid = sum2_bid/numpnts(bid_wave)
	variable std_bid = sqrt(mean2_bid - mean1_bid*mean1_bid)
	mean1_ask = sum1_ask/numpnts(ask_wave)
	mean2_ask = sum2_ask/numpnts(ask_wave)
	variable std_ask = sqrt(mean2_ask - mean1_ask*mean1_ask)
	variable sum_co = 0
	for(i=0;i<numpnts(bid_wave);i+=1)
		sum_co += (bid_wave[i] - mean1_bid) * (ask_wave[i] - mean1_ask)
	endfor
	variable result = (sum_co/numpnts(bid_wave))/(std_bid * std_ask)
	print("The Correlation Between waves is " + num2str(result))
End

Function Correlation()
	DFREF now_dfr = getDataFolderDFR()
	DFREF working_dfr = get_Correlation_datafolder()
	string Window_name = "Correlation_Controler"
	string DT_BID_BUTTON_NAME = "DT_Bid_Button"
	string DT_ASK_BUTTON_NAME = "DT_Ask_Button"
	if(wintype(Window_name) != 0)
		Dowindow /F $(Window_name)
	else
		setdatafolder working_dfr			// Make global variables on working directory
		string /G DT_Global_WaveName_Bid = ""
		string /G DT_Global_WaveName_Ask = ""
		setdatafolder now_dfr
		NewPanel /W=(50,50,240,280) /K=1 as Window_name
		RenameWindow $(S_name), $(Window_name)
		GroupBox index_Group1, pos={10,5},size={170,50},title="select wave"
		Button $(DT_BID_BUTTON_NAME), pos={20,25},size={150,20}
		MakeButtonIntoWSPopupButton(Window_name,DT_BID_BUTTON_NAME,"CL_PopupWaveSelector_Bid",options=PopupWS_OptionFloat)
		GroupBox index_Group2, pos={10,55},size={170,50},title="select wave"
		Button $(DT_ASK_BUTTON_NAME),pos={20,75},size={150,20}
		MakeButtonIntoWSPopupButton(Window_name,DT_ASK_BUTTON_NAME,"CL_PopupWaveSelector_Ask",options=PopupWS_OptionFloat)
		//checkbox DT_checkbox_only_Bid pos={10,100},title="If you wanna use bid, please check"
		
		Button DT_execute_Button, pos={30,120},size={130,20},proc=Correlation_core,title="Execute", fSize=10
		
		GroupBox index_group3, pos={10,150},size={170,0},title="Correlation Value is displayed"
		GroupBox index_group4,pos={10,170},size={170,0},title=" on Terminal"
	endif
End


// ----------- H2 room analize -------------
Function three_points_positive_peak(wave_name)
	string wave_name
	string output_name = "index_of_peaks"
	make /O /N=(0) $(output_name) /wave=out
	wave w = $(wave_name)
	variable i
	variable a,b,c
	//variable height = 6
	for(i=0;i<numpnts(w)-2;)
		a = w[i]
		b = w[i+1]
		c = w[i+2]
		if(a < b && b > c)
			if(b > 9)
				make /O /N=(numpnts(out)+1) $(output_name)
				out[numpnts(out)-1] = i+1
			endif
			i += 2
		elseif(a == b && b > c)
			if(b > 9)
				make /O /N=(numpnts(out)+1) $(output_name)
				out[numpnts(out)-1] = i+1
			endif
			i += 2
		elseif(a > b && b > c)
			i += 2
		else
			i += 1
		endif
	endfor
End

Function three_points_negative_peak(wave_name)
	string wave_name
	string output_name = "index_of_peaks_negaive"
	make /O /N=(0) $(output_name) /wave=out
	wave w = $(wave_name)
	variable i
	variable a,b,c
	variable height = 3
	for(i=0;i<numpnts(w)-2;)
		a = w[i]
		b = w[i+1]
		c = w[i+2]
		if(a > b && b < c)
			if(b < 1 && c - b > 4)
				make /O /N=(numpnts(out)+1) $(output_name)
				out[numpnts(out)-1] = i+1
			endif
			i += 2
		elseif(a == b && b < c)
			if(b < 1 && c - b > 4)
				make /O /N=(numpnts(out)+1) $(output_name)
				out[numpnts(out)-1] = i+1
			endif
			i += 2
		elseif(a < b && b < c)
			i += 2
		else
			i += 1
		endif
	endfor
End


#if (IgorVersion() < 7) 

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
		NewPanel/W = (50,50,240,700)/K=1 as "AutoFit_Panel"
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
		checkbox NS_autoFit_checkbox pos={10,82},title="fix on / off (for Freq)"
		//checkbox NS_autoFit_checkbox_area pos={10,102},title = "fix on / off (for Area)"
		SetVariable setwave_name pos={10,122},size={150,20},value=AF_num_of_peaks,format="%d",title="num of peaks :"
		string tmp_folder = "\"" + tmp_str + "\""
		PopupMenu NF_select_folder,mode=1,Value=#tmp_folder,title="select waves folder",pos={10,146},size={17,20},proc=PopupMenuAction_for_autofit
		Button draw_button_NS_AutoFit, pos={10,176},size={170,20},proc=NS_Auto_Fit_button,title="Do Fit", fSize=10
		
		variable space = 20
		Groupbox option_group1,pos={10,200},size={170,220},title="select parameters to hold (area)"
		checkbox NS_autoFit_checkbox_area_1 pos={10,200+space},title="1 gauss for area"
		checkbox NS_autoFit_checkbox_area_2 pos={10,200+2*space},title="2 gauss for area"
		checkbox NS_autoFit_checkbox_area_3 pos={10,200+3*space},title="3 gauss for area"
		checkbox NS_autoFit_checkbox_area_4 pos={10,200+4*space},title="4 gauss for area"
		checkbox NS_autoFit_checkbox_area_5 pos={10,200+5*space},title="5 gauss for area"
		checkbox NS_autoFit_checkbox_area_6 pos={10,200+6*space},title="6 gauss for area"
		checkbox NS_autoFit_checkbox_area_7 pos={10,200+7*space},title="7 gauss for area"
		checkbox NS_autoFit_checkbox_area_8 pos={10,200+8*space},title="8 gauss for area"
		checkbox NS_autoFit_checkbox_area_9 pos={10,200+9*space},title="9 gauss for area"
		checkbox NS_autoFit_checkbox_area_10 pos={10,200+10*space},title="10 gauss for area"
		
		Groupbox option_group2,pos={10,420},size={170,220},title="---------------------- (width)"
		checkbox NS_autoFit_checkbox_width_1 pos={10,420+space},title="1 gauss for width"
		checkbox NS_autoFit_checkbox_width_2 pos={10,420+2*space},title="2 gauss for width"
		checkbox NS_autoFit_checkbox_width_3 pos={10,420+3*space},title="3 gauss for width"
		checkbox NS_autoFit_checkbox_width_4 pos={10,420+4*space},title="4 gauss for width"
		checkbox NS_autoFit_checkbox_width_5 pos={10,420+5*space},title="5 gauss for width"
		checkbox NS_autoFit_checkbox_width_6 pos={10,420+6*space},title="6 gauss for width"
		checkbox NS_autoFit_checkbox_width_7 pos={10,420+7*space},title="7 gauss for width"
		checkbox NS_autoFit_checkbox_width_8 pos={10,420+8*space},title="8 gauss for width"
		checkbox NS_autoFit_checkbox_width_9 pos={10,420+9*space},title="9 gauss for width"
		checkbox NS_autoFit_checkbox_width_10 pos={10,420+10*space},title="10 gauss for width"
		
		
		
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
	
	//switch(i)
	//	case 4:
	//		return "0010"
	//	case 5:
	//		return "00100"
	//	case 8:
	//		return "00100010"
	//	case 11:
	//		return "00100010010"
	//	case 14:
	//		return "00100010010010"
	//	case 17:
	//		return "00100010010010010"
	//	case 20:
	//		return "00100010010010010010"
	//	case 23:
	//		return "00100010010010010010010"
	//	case 26:
	//		return "00100010010010010010010010"
	//	case 29:
	//		return "00100010010010010010010010010"
	//	case 32:
	//		return "00100010010010010010010010010010"
	//	default:
	//		doalert 0,"error of Fitting You should reduce waves"
	//		return ""
	//endswitch
	return NS_GET_HOLD_STR(i)
End

static Function /s NS_GET_HOLD_STR(i)		//Freq fix したときのみなのであとでarea と widthだけ固定できるよう変更
	variable i
	
	variable area_1
	variable area_2
	variable area_3
	variable area_4
	variable area_5
	variable area_6
	variable area_7
	variable area_8
	variable area_9
	variable area_10
	
	variable width_1
	variable width_2
	variable width_3
	variable width_4
	variable width_5
	variable width_6
	variable width_7
	variable width_8
	variable width_9
	variable width_10
	
	
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_area_1
	area_1 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_area_2
	area_2 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_area_3
	area_3 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_area_4
	area_4 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_area_5
	area_5 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_area_6
	area_6 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_area_7
	area_7 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_area_8
	area_8 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_area_9
	area_9 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_area_10
	area_10 = V_Value
	
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_width_1
	width_1 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_width_2
	width_2 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_width_3
	width_3 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_width_4
	width_4 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_width_5
	width_5 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_width_6
	width_6 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_width_7
	width_7 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_width_8
	width_8 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_width_9
	width_9 = V_Value
	controlinfo /W = $("AutoFit_Panel_p") NS_autoFit_checkbox_width_10
	width_10 = V_Value
	
	switch(i)
		case 4:
			return "0" + num2str(area_1) + "1" + num2str(width_1)
		case 5:
			return "0" + num2str(area_1) + "1" + num2str(width_1) + "0"
		case 8:
			return "0" + num2str(area_1) + "1" + num2str(width_1) + "0" + num2str(area_2) + "1" + num2str(width_2)
		case 11:
			return "0" + num2str(area_1) + "1" + num2str(width_1) + "0" + num2str(area_2) + "1" + num2str(width_2) + num2str(area_3) + "1" + num2str(width_3) 
		case 14:
			return 	"0" + num2str(area_1) + "1" + num2str(width_1) + "0" + num2str(area_2) + "1" + num2str(width_2) + num2str(area_3) + "1" + num2str(width_3) + num2str(area_4) + "1" + num2str(width_4)
		case 17:
			string tt = 	"0" + num2str(area_1) + "1" + num2str(width_1) + "0" + num2str(area_2) + "1" + num2str(width_2) + num2str(area_3) + "1" + num2str(width_3) + num2str(area_4) + "1" + num2str(width_4)
			tt += num2str(area_5) + "1" + num2str(width_5)
			return tt
		case 20:
			string tt2 = 	"0" + num2str(area_1) + "1" + num2str(width_1) + "0" + num2str(area_2) + "1" + num2str(width_2) + num2str(area_3) + "1" + num2str(width_3) + num2str(area_4) + "1" + num2str(width_4)
			tt2 += num2str(area_5) + "1" + num2str(width_5) + num2str(area_6) + "1" + num2str(width_6)
			return tt2
		case 23:
			string tt3 = 	"0" + num2str(area_1) + "1" + num2str(width_1) + "0" + num2str(area_2) + "1" + num2str(width_2) + num2str(area_3) + "1" + num2str(width_3) + num2str(area_4) + "1" + num2str(width_4)
			tt3 += num2str(area_5) + "1" + num2str(width_5) + num2str(area_6) + "1" + num2str(width_6) + num2str(area_7) + "1" + num2str(width_7)
			return tt3
		case 26:
			string tt4 = 	"0" + num2str(area_1) + "1" + num2str(width_1) + "0" + num2str(area_2) + "1" + num2str(width_2) + num2str(area_3) + "1" + num2str(width_3) + num2str(area_4) + "1" + num2str(width_4)
			tt4 += num2str(area_5) + "1" + num2str(width_5) + num2str(area_6) + "1" + num2str(width_6) + num2str(area_7) + "1" + num2str(width_7) + num2str(area_8) + "1" + num2str(width_8)
			return tt4
		case 29:
			string tt5 = 	"0" + num2str(area_1) + "1" + num2str(width_1) + "0" + num2str(area_2) + "1" + num2str(width_2) + num2str(area_3) + "1" + num2str(width_3) + num2str(area_4) + "1" + num2str(width_4)
			tt5 += num2str(area_5) + "1" + num2str(width_5) + num2str(area_6) + "1" + num2str(width_6) + num2str(area_7) + "1" + num2str(width_7) + num2str(area_8) + "1" + num2str(width_8) + num2str(area_9) + "1" + num2str(width_9)
			return tt5
		case 32:
			string tt6 = 	"0" + num2str(area_1) + "1" + num2str(width_1) + "0" + num2str(area_2) + "1" + num2str(width_2) + num2str(area_3) + "1" + num2str(width_3) + num2str(area_4) + "1" + num2str(width_4)
			tt6 += num2str(area_5) + "1" + num2str(width_5) + num2str(area_6) + "1" + num2str(width_6) + num2str(area_7) + "1" + num2str(width_7) + num2str(area_8) + "1" + num2str(width_8) + num2str(area_9) + "1" + num2str(width_9)
			tt6 += num2str(area_10) + "1" + num2str(width_10)
			return tt6
		
		
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

static Function autoFit(x0,x1,W_coef,dfr,out_dfr,length,wave_names,[rev,fix,subtract,fix2]) //length is the number of peaks
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


//--------------------------------------- New Hoitaku's AutoFit  --------------------------------

Static Function BinarySearchClipped(w,x)
	WAVE w
	Variable x
	
	Variable p= BinarySearch(w,x)
	if( p == -2 )
		p= numpnts(w)-1
	elseif( p == -1 )
		p= 0
	endif
	
	return p
End

static Function MPF2_SetDataPointRange(gname, YData, XData, RangeBegin, RangeEnd, RangeReversed)
	string gname
	Wave YData
	Wave/Z XData
	Variable &RangeBegin
	Variable &RangeEnd
	Variable &RangeReversed

	RangeBegin= 0;
	RangeEnd= numpnts(YData)-1;
	
	Variable te= RangeEnd
	Variable V_Flag= 0

	CheckDisplayed/W=$gname YData
	Variable isGraphed= V_Flag
	Variable checkAxis = 0
	if( isGraphed )
		ControlInfo/W=$(gname+"#MultiPeak2Panel") MPF2_UserCursorsCheckbox
		if (1) //always on
			if (strlen(CsrInfo(A, gname)) == 0)
				DoAlert 0, "The Use Graph Cursors checkbox is checked, but the A cursor is not on the graph."
				checkAxis = 1
			endif
			if (strlen(CsrInfo(B, gname)) == 0)
				DoAlert 0, "The Use Graph Cursors checkbox is checked, but the B cursor is not on the graph."
				checkAxis = 1
			endif
			if (checkAxis == 0)
				RangeBegin= pcsr(A, gname)
				RangeEnd= pcsr(B, gname)
			endif
		else
			checkAxis = 1
		endif
		if (checkAxis)
			GetAxis /Q bottom
			if(!WaveExists(XData))
				RangeBegin= max(x2pnt(YData,V_min), 0)
				RangeEnd= min(x2pnt(YData,V_max), numpnts(YData)-1)
			else
				RangeBegin=BinarySearchClipped(XData,V_min)
				RangeEnd=BinarySearchClipped(XData,V_max)
			endif
		endif
	endif
	RangeReversed= RangeBegin>RangeEnd
	if( RangeReversed )
		variable tmp= RangeBegin
		RangeBegin= RangeEnd
		RangeEnd= tmp
	endif
	return 0
End

static Function/S MPF2_FolderNameFromSetNumber(setnumber)
	Variable setnumber
	
	return "MPF_SetFolder_"+num2str(setnumber)
end

static Function/S MPF2_FolderPathFromSetNumber(setnumber)
	Variable setnumber
	
	return "root:Packages:MultiPeakFit2:"+MPF2_FolderNameFromSetNumber(setnumber)
end

Static Function GetSetNumberFromWinName(windowName)
	String windowName
	
	String windowWithData
	
	Variable poundPos = strsearch(windowName, "#", 0)
	if (poundPos < 0)
		windowWithData = windowName
	else
		poundPos = strsearch(windowName, "#", poundPos+1)
		if (poundPos < 0)
			windowWithData = windowName
		else
			windowWithData = windowName[0,poundPos-1]
		endif
	endif
	
	return str2num(GetUserData(windowWithData, "", "MPF2_DataSetNumber"))
end

Static Function/S MPF2_HoldStringForPeakListItem(theItem, DatafolderPath, thePanelWin)
	String theItem, DatafolderPath, thePanelWin
	
	String SaveDF = GetDataFolder(1)
	SetDataFolder DatafolderPath
	
	Wave/T HoldStrings
	if (!WaveExists(HoldStrings))
		SetDataFolder SaveDF
		return ""
	endif
	
	String holdString = ""
	Variable isBaseLine = CmpStr(theItem, "Baseline") == 0
	String children=""
	Variable i
	Variable numItems
	Variable rownumber
	Variable HoldStringsRow
	
	// It is possible to get into this function without having set the size of the hold strings wave correctly
	Wave/Z wpi = W_AutoPeakInfo
	Variable npeaks = 0
	if (WaveExists(wpi))
		npeaks = DimSize(wpi, 0)
	endif
	if (DimSize(HoldStrings, 0) < npeaks+1)
		Variable oldSize = DimSize(HoldStrings, 0)
		Redimension/N=(npeaks+1) HoldStrings
		HoldStrings[oldSize, npeaks] = ""
	endif
	
	if (isBaseLine)
		HoldStringsRow = 0
	else	
		sscanf theItem, "Peak %d", HoldStringsRow
		HoldStringsRow += 1							// to account for the fact that the baseline holds are in row 0
	endif

	rownumber = WMHL_GetRowNumberForItem(thePanelWin, "MPF2_PeakList", theItem)
	if (WMHL_RowIsOpen(thePanelWin, "MPF2_PeakList", rownumber))
		children = WMHL_ListChildRows(thePanelWin, "MPF2_PeakList", rownumber)
		numItems = ItemsInList(children)
		for (i = 0; i < numitems; i += 1)
			rownumber = str2num(StringFromList(i, children))
			Variable SelWaveValue = WMHL_GetExtraColumnSelValue(thePanelWin, "MPF2_PeakList", 1, rownumber)
			if (SelWaveValue & 0x10)
				holdString += "1"
			else
				holdString += "0"
			endif
		endfor
	else
		holdString = HoldStrings[HoldStringsRow]
	endif
	
	Variable numchars = strlen(holdString)
	Variable Char1 = char2num("1")
	
	for (i = 0; i < numchars; i += 1)
		if (char2num(holdString[i]) == Char1)
			break;			// found a 1
		endif
	endfor
	if (i == numchars)
		holdString = ""		// didn't find a 1
	endif
	
	SetDataFolder SaveDF
	return holdString
end


Static Function MPF2_RefreshHoldStrings(PanelWin)
	String PanelWin
	
	Variable setNumber = GetSetNumberFromWinName(PanelWin)
	String DFpath = MPF2_FolderPathFromSetNumber(setNumber)
	Wave/T HoldStrings = $(DFPath+":HoldStrings")
	String listPanel = PanelWin+"#P1"

	Variable rownumber = WMHL_GetRowNumberForItem(listPanel, "MPF2_PeakList", "Baseline")
	if (WMHL_RowIsOpen(listPanel, "MPF2_PeakList", rownumber))
		HoldStrings[0] = MPF2_HoldStringForPeakListItem("Baseline", DFPath, listPanel)
	endif

	Variable i=0
	do
		String peakItem = "Peak "+num2str(i)
		rownumber = WMHL_GetRowNumberForItem(listPanel, "MPF2_PeakList", peakItem)
		if (rownumber < 0)
			break;
		endif
		HoldStrings[i+1] = MPF2_HoldStringForPeakListItem(peakItem, DFPath, listPanel)
		
		i += 1
	while(1)
end

Static Function CreateCWavesInCDFFromAutoPkInfo(AutoPeakInfo, peakTypeName, coefWaveNameFormat)
	Wave AutoPeakInfo
	String peakTypeName
	String coefWaveNameFormat
	
	Variable npeaks = DimSize(AutoPeakInfo, 0)
	
	FUNCREF MPF2_FuncInfoTemplate infoFunc=$(peakTypeName+PEAK_INFO_SUFFIX)
	Variable nparams
	String GaussGuessConversionFuncName = infoFunc(PeakFuncInfo_GaussConvFName)
	if (strlen(GaussGuessConversionFuncName) == 0)
	else
		FUNCREF MPF2_GaussGuessConvTemplate gconvFunc=$GaussGuessConversionFuncName
	endif
	
	String newWName
	
	Variable i
	for (i = 0; i < npeaks; i += 1)
		sprintf newWName, coefWaveNameFormat, i
		Make/D/O/N=(DimSize(AutoPeakInfo, 1)) $newWName
		Wave w = $newWName
		w = AutoPeakInfo[i][p]
		gconvFunc(w)
	endfor
end

Static Function CreateCoefWavesFromAutoPeakInfo(setNumber, AutoPeakInfo, peakTypeName)
	Variable setNumber
	Wave AutoPeakInfo
	String peakTypeName
	
//	NVAR currentSetNumber = root:Packages:MultiPeakFit2:currentSetNumber
	String DFpath = MPF2_FolderPathFromSetNumber(setNumber)
	String saveDF = GetDataFolder(1)
	SetDataFolder DFPath
	
	CreateCWavesInCDFFromAutoPkInfo(AutoPeakInfo, peakTypeName, "Peak %d Coefs")
	
	SetDataFolder saveDF
end

Static Function MPF2_PutAutoPeakResultIntoList(setNumber, autoPeakInfo, initializeBaseline [, listOfPeakTypes])
	Variable setNumber
	Wave autoPeakInfo
	Variable initializeBaseline
	String listOfPeakTypes
	
	String DFpath = MPF2_FolderPathFromSetNumber(setNumber)
	SVAR gname = $(DFpath+":GraphName")
	SVAR YWvName = $(DFpath+":YWvName")
	SVAR XWvName = $(DFpath+":XWvName")
	Wave YData = $YWvName
	Wave/Z XData = $XWvName
	
	Variable i, theRow
	if (!WaveExists(autoPeakInfo) || (DimSize(autoPeakInfo, 0) == 0))
		// If there's any peaks in the listbox get rid of them
		if (!initializeBaseline)
			i = 0
			do
				theRow = WMHL_GetRowNumberForItem(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Peak "+num2str(i))
				if (theRow < 0)
					break;
				endif
				WMHL_DeleteRowAndChildren(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", theRow)
				i += 1
			while (1)
		endif
		
		return 0
	endif
	
	Variable nPeaks = DimSize(	autoPeakInfo, 0)
	Variable currentRow = 0
	Variable reOpenBaseline = 0
	
	if (initializeBaseline)
		// Start from scratch
		ControlInfo/W=$gname#MultiPeak2Panel#P1 MPF2_PeakList
		Variable listHeight = V_height
		Variable listWidth = V_width
		Variable listTop = V_top
		Variable listLeft = V_left
		
		Variable BaselineRow = WMHL_GetRowNumberForItem(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Baseline")
		String baselineStr
		if (BaselineRow == 0)
			baselineStr = WMHL_GetExtraColumnData(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0, BaselineRow)
		else
			baselineStr = "Constant"+MENU_ARROW_STRING
		endif

		ListBox MPF2_PeakList win=$gname#MultiPeak2Panel#P1,userdata(MPF2_DataSetNumber)=num2str(setnumber)
		MakeListIntoHierarchicalList(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "PeakListOpenNotify", selectionMode=WMHL_SelectionContinguous, userListProc="MPF2_PeakListProc")
		WMHL_SetNotificationProc(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "PeakListClosingNotify", WMHL_SetClosingNotificationProc)
		WMHL_AddColumns(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", 6)
		
		WMHL_AddObject(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "", "Baseline", 1)
		WMHL_ExtraColumnData(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0, currentRow, baselineStr, 0)
		currentRow += 1
	else
		// preserve the baseline info. Just delete all the peaks and add them back again.
		reOpenBaseline = WMHL_RowIsOpen(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0)
		WMHL_CloseAContainer(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Baseline")
		i = 0
		do
			theRow = WMHL_GetRowNumberForItem(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Peak "+num2str(i))
			if (theRow < 0)
				break;
			endif
			WMHL_DeleteRowAndChildren(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", theRow)
			i += 1
		while (1)
		currentRow = 1
	endif
	
	String peakTypeName = "Gauss"
	for (i = 0; i < npeaks; i += 1)
		WMHL_AddObject(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "", "Peak "+num2str(i), 1)
		if (!ParamIsDefault(listOfPeakTypes))
			peakTypeName = StringFromList(i+1, listOfPeakTypes)		// i+1 because the first is the baseline type
		endif
		WMHL_ExtraColumnData(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0, currentRow, peakTypeName+MENU_ARROW_STRING, 0)
		currentRow += 1
	endfor
	
	ControlInfo/W=gname#MultiPeak2Panel#P1 MPF2_PeakList
	
	Variable updateRejected = NumVarOrDefault("root:Packages:MultiPeakFit2:updateRejected", 0)
	if (!updateRejected)
		NVAR MPF2ConstraintsShowing = $(DFPath+":MPF2ConstraintsShowing")
		if (MPF2ConstraintsShowing)
			ListBox MPF2_PeakList win=$gname#MultiPeak2Panel#P1,widths={4,12,15,11,5,10,6,10}
		else 
			ListBox MPF2_PeakList win=$gname#MultiPeak2Panel#P1,widths={4,12,15,11,0,0,0,0}
		endif
	else	
		ListBox MPF2_PeakList win=$gname#MultiPeak2Panel#P1,widths={4,15,16,10}
	endif
	
	if (reOpenBaseline)
		WMHL_OpenAContainer(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Baseline")
	endif
end

static Function MPF2_getFitCurvePoints(panelName)
	String panelName
	
	Variable points = 0
	
	ControlInfo/W=$(panelName+"#P3") MPF2_SetFitCurvePoints
	String pointsStr = S_value
	if (CmpStr(S_value, "Auto") != 0)
		points = str2num(pointsStr)
		if (Numtype(points) != 0)
			points = 0
		endif
	endif
	
	return points
end
 
 static Function HO_MPF2_PeakResultsButtonProc(graphStr,setNumber)
	string graphStr
	variable setNumber
	
//	String gname = WinName(0,1)
	MPF2_DoPeakResults(GetSetNumberFromWinName(graphStr))
	String MPF2_ResultsPanelName = "MPF2_ResultsPanel"+"_"+num2str(setNumber)
	Dowindow /K  $(MPF2_ResultsPanelName)
	//Killwindow /Z $("Multipeak Fit 2 Results")
End

 
 // main procedure
 
static Function HO_MPF2_DoFitButtonProc(setNumber,windowname)
	variable setNumber
	string windowname
	
	STRUCT MPFitInfoStruct MPStruct
	String DFpath = MPF2_FolderPathFromSetNumber(setNumber)
	
	Wave/Z wpi = $(DFPath+":"+"W_AutoPeakInfo")
	if (!WaveExists(wpi))
		DoAlert 0, "There are no peaks to fit."
		return -1
	endif
	MPStruct.NPeaks = DimSize(wpi, 0)
	SVAR YWvName = $(DFpath+":YWvName")
	SVAR XWvName = $(DFpath+":XWvName")
	Wave MPStruct.yWave = $YWvName
	Wave/Z MPStruct.xWave = $XWvName
	SVAR gname = $(DFpath+":GraphName")
	String listPanelName = gname+"#MultiPeak2Panel#P1"
	
	Wave/Z MPStruct.weightWave = $PopupWS_GetSelectionFullPath( gname+"#MultiPeak2Panel#P3", "MPF2_SelectWeightWave")
	Wave/Z MPStruct.maskWave = $PopupWS_GetSelectionFullPath( gname+"#MultiPeak2Panel#P3", "MPF2_SelectMaskWave")

	String saveDF = GetDataFolder(1)
	SetDataFolder DFPath

	//// do a check on the Use Graph Cursors checkbox to determine xPointRangeBegin&End
	String multiPeakPanelName = gname+"#MultiPeak2Panel"
	if (1)			// use cursors
		NVAR XPointRangeBegin
		NVAR XPointRangeEnd	
		MPStruct.XPointRangeBegin = XPointRangeBegin
		MPStruct.XPointRangeEnd = XPointRangeEnd
	else					// don't use cursors - use visible graph range size
		Variable RangeBegin, RangeEnd, RangeReversed
		MPF2_SetDataPointRange(gname, MPStruct.yWave, MPStruct.xWave, RangeBegin, RangeEnd, RangeReversed)
		MPStruct.XPointRangeBegin = RangeBegin
		MPStruct.XPointRangeEnd = RangeEnd
	endif

	NVAR MPF2_FitCurvePoints
	MPStruct.FitCurvePoints = MPF2_FitCurvePoints

	String/G FuncListString=""
	
	Variable nBLParams
	String ParamNameList
	String pwname
	Variable i
	
	String OneHoldString = ""
	Wave/T HoldStrings
	
	Variable BaselineRow = WMHL_GetRowNumberForItem(listPanelName, "MPF2_PeakList", "Baseline")
	String baselineStr = WMHL_GetExtraColumnData(listPanelName, "MPF2_PeakList", 0, BaselineRow)
	MPStruct.ListOfFunctions = MPF2_PeakOrBLTypeFromListString(baselineStr)+";"
	//MPStruct.ListOfFunctions = MPF2_PeakOrBLTypeFromListString("Linear\JR\W523")+";"
	Variable doBaseLine = CmpStr(baselineStr, "None"+MENU_ARROW_STRING) != 0
	MPStruct.ListOfCWaveNames = "Baseline Coefs;"		// if baseline type is "None", this wave probably doesn't exist, but it doesn't matter because it will be ignored
	if (doBaseLine)
		MPStruct.ListOfHoldStrings = MPF2_HoldStringForPeakListItem("Baseline", DFpath, listPanelName)+";"
	else
		MPStruct.ListOfHoldStrings = ";"
	endif
	
	for (i = 0; i < MPStruct.NPeaks; i += 1)
		MPStruct.ListOfCWaveNames += "Peak "+num2istr(i)+" Coefs;"
		String peakItem = "Peak "+num2istr(i)
		Variable theRow = WMHL_GetRowNumberForItem(listPanelName, "MPF2_PeakList", peakItem)
		String PeakTypeName = MPF2_PeakOrBLTypeFromListString( WMHL_GetExtraColumnData(listPanelName, "MPF2_PeakList", 0, theRow) )
		MPStruct.ListOfFunctions += PeakTypeName+";"
		MPStruct.ListOfHoldStrings += MPF2_HoldStringForPeakListItem(peakItem, DFpath, listPanelName)+";"
	endfor
//print "Function list = ["+FuncListString+"]"
//print "Function list has ", strlen(FuncListString), "characters"

	//// Check inter-peak constraints ////
	NVAR/Z MPF2ConstraintsShowing = $(DFPath+":MPF2ConstraintsShowing")
	if (!NVAR_Exists(MPF2ConstraintsShowing))
		Variable/G $(DFPath+":MPF2ConstraintsShowing")
		NVAR MPF2ConstraintsShowing = $(DFPath+":MPF2ConstraintsShowing")
		MPF2ConstraintsShowing = 0
	endif
	if (MPF2ConstraintsShowing)	
		SVAR /Z interPeakString = $(DFPath+":interPeakConstraints")
		if (!SVAR_exists(interPeakString))
			String /G $(DFPath+":interPeakConstraints")
			SVAR interPeakString = $(DFPath+":interPeakConstraints")
		endif

		String notebookName = StringFromList(0,windowname,"#")+"#"+StringFromList(1,windowname,"#")+"#P3#MPF2_InterPeakConstraints"
		Notebook $notebookName, selection={startOfFile, endOfFile}
		GetSelection notebook, $notebookName, 2
		interPeakString = S_Selection[0, strlen(S_Selection)-1]
		interPeakString = ReplaceString("\n", interPeakString, "")
		interPeakString = ReplaceString("\r", interPeakString, "")
		Notebook $notebookName, visible=0
		Notebook $notebookName, visible=1
			
		Variable isValid = MPF2_ValidateConstraint(setNumber, windowname, interPeakString)
		if (!isValid)
			return -1
		endif

		Wave /T MPStruct.constraints = getGlobalConstraintsWave(setNumber)

		Duplicate /O MPStruct.constraints, $(DFPath+":MPF2_ConstraintsBackup")
	else	
		KillWaves /Z $(DFPath+":MPF2_ConstraintsBackup")
	endif

	// Added to have a FuncFit ready constraints wave ready.  

	MPStruct.fitOptions = 4
	
	MPF2_SaveFunctionTypes(gname+"#MultiPeak2Panel")
	
//Variable etime = ticks	
	MPF2_DoMPFit(MPStruct, DFPath+":")
//etime = ticks-etime
//print "Time for fit: ", etime/60," seconds"	
	FuncListString = MPStruct.FuncListString

	String alertMsg
//	if (MPStruct.fitError || ((MPStruct.fitQuitReason > 0) && (MPStruct.fitQuitReason < 3)) )
//		Variable doRestore = 1
//		if (MPStruct.fitQuitReason == 2)
//			alertMsg = "Multi-peak Fit cancelled."
//		else
//			alertMsg = "Multi-peak  fit failed:"
//		endif
//		if (MPStruct.fitError & 2)
//			alertMsg += "\rSingular matrix error."
//		endif
//		if (MPStruct.fitError & 4)
//			alertMsg += "\rOut of memory."
//		endif
//		if (MPStruct.fitError & 8)
//			alertMsg += "\rFunction return NaN or INF."
//		endif
//		if (MPStruct.fitQuitReason == 1)
//			alertMsg += "\rNumber of iterations exceded the limit."
//			doRestore = 0		// allow the fit to continue from where it left off if Do Fit is clicked again
//		endif
//		if (doRestore)
//			MPF2_RestoreCoefWavesFromBackup(MPStruct.ListOfCWaveNames, DFPath)
//		endif
//		DoAlert 0, alertMsg
	Variable doRestore = 0
	if (MPStruct.fitError || MPstruct.fitQuitReason)
		doRestore = 1
		if (MPStruct.fitError)
			DoAlert 0, "Multi-peak Fit failed: \r\r"+MPstruct.fitErrorMsg
		else
			switch (MPstruct.fitQuitReason)
				case 1:
					DoAlert 0, "Multi-peak fit exceded the iteration limit. Click Fit again to continue."
					doRestore = 0
					break;
				case 2:
					DoAlert 0, "Multi-peak fit cancelled."
					break;
				case 3:
					DoAlert 0, "Multi-peak fit not progressing. Chances are the fit is good."
					doRestore = 0
					break;
			endswitch
		endif	
		if (doRestore)
			MPF2_RestoreCoefWavesFromBackup(MPStruct.ListOfCWaveNames, DFPath)
		endif
	endif
	if (doRestore == 0)
		Variable/G MPF2_FitDate = MPStruct.dateTimeOfFit
		Variable/G MPF2_FitPoints = MPStruct.fitPnts
		Variable/G MPF2_FitChiSq = MPStruct.chisq
		// Now update the list with the fit results
		if (doBaseLine)
			if (WMHL_RowIsOpen(listPanelName, "MPF2_PeakList", BaselineRow))
				Wave 'Baseline Coefs'
				String baselineRows = WMHL_ListChildRows(listPanelName, "MPF2_PeakList", BaselineRow)
				Variable numBLRows = ItemsInList(baseLineRows)
				for (i = 0; i < numBLRows; i += 1)
					Variable BLcoefRow = str2num(stringFromList(i, baselineRows))
					WMHL_ExtraColumnData(listPanelName, "MPF2_PeakList", 0, BLcoefRow, num2str('Baseline Coefs'[i]), 1)
				endfor
			endif
		endif

		for (i = 0; i < MPStruct.NPeaks; i += 1)
			Variable PeakRow = WMHL_GetRowNumberForItem(listPanelName, "MPF2_PeakList", "Peak "+num2istr(i))
			if (WMHL_RowIsOpen(listPanelName, "MPF2_PeakList", PeakRow))
				Wave coefs = $("Peak "+num2istr(i)+" Coefs")
				String coefRows = WMHL_ListChildRows(listPanelName, "MPF2_PeakList", PeakRow)
				Variable numCoefRows = ItemsInList(coefRows)
				Variable j
				for (j = 0; j < numCoefRows; j += 1)
					Variable peakCoefRow = str2num(stringFromList(j, coefRows))
					WMHL_ExtraColumnData(listPanelName, "MPF2_PeakList", 0, peakCoefRow, num2str(coefs[j]), 1)
				endfor
			endif
//			Make/N=5/O tempAutoPickInfo
//			MPF2_GetSimulatedAutoPickData(PeakRow, listPanelName, tempAutoPickInfo)
//			wpi[i][] = tempAutoPickInfo[q]
		endfor
		
		MPF2_RefreshPeakResults(setNumber)
	endif

	NVAR negativePeaks = negativePeaks
	NVAR displayPeaksFullWidth = $(DFpath+":displayPeaksFullWidth")
	MPF2_AddPeaksToGraph(setNumber, wpi, 1, 1, displayPeaksFullWidth)
	MPF2_AddFitCurveToGraph(setNumber, wpi, MPStruct.yWave, MPStruct.xWave, 1, overridePoints=MPF2_getFitCurvePoints(gname+"#MultiPeak2Panel"))
	
	SetDataFolder saveDF		
End


//この関数で一つ一つFitting Operation
static Function NAF_Element_handler(backup_wd,working_dfr,NAF_DFR,waves_dfr,y_wave_name,x_wave,hold_wave,guessing_wave,num_of_gaussian,Acursor,Bcursor)
	DFREF backup_wd,working_dfr
	DFREF NAF_DFR,waves_dfr
	string y_wave_name
	wave x_wave
	wave hold_wave,guessing_wave
	variable num_of_gaussian,Acursor,Bcursor
	
	variable i
	
	DFREF local_DFR = getDataFolderDFR()
	setdatafolder NAF_DFR
		SVar N_AF_Baseline_Func_Name
		Svar N_AF_Fitting_Func_Name
	setdatafolder local_DFR
	
	wave /Z /SDFR=waves_dfr y_wave = $(y_wave_name)
	if(waveexists(y_wave) == 0)
		setdatafolder working_dfr
		return -1
	endif
	
	MPF2_StartNewMPFit(0,"New Graph",GetWavesDataFolder(y_wave,2),GetWavesDataFolder(x_wave,2),1,2)
	Nvar currentSetNumber = root:Packages:MultiPeakFit2:currentSetNumber
	string graphStr = WinName(0,1)
	DFREF current_working_dfr = $(MPF2_FolderPathFromSetNumber(currentSetNumber))
	cursor /A=1 /W=$(graphStr) A $(nameofwave(y_wave)) Acursor
	cursor /A=1 /W=$(graphStr) B $(nameofwave(y_wave)) Bcursor
	setDataFolder current_working_dfr
	MPF2_SetDataPointRange(graphStr,y_wave,x_wave,Acursor,Bcursor,Acursor)
	Variable/G XPointRangeBegin = Acursor
	Variable/G XPointRangeEnd = Bcursor
	Variable/G XPointRangeReversed = 0
	
	//SVAR YWvName = $(DFpath+":YWvName")
	//SVAR XWvName = $(DFpath+":XWvName")
	//Wave yw = $YWvName
	//Wave/Z xw = $XWvName
	
	Make /D/N=(num_of_gaussian,5)/O $("W_AutoPeakInfo") /wave=wpi
	
	setdatafolder NAF_DFR
	Svar N_AF_Baseline_Func_Name
	setdatafolder working_dfr
	
	if(cmpstr("Constant",N_AF_Baseline_Func_Name) == 0)
		//baseline const
		wave tmp = $("Baseline Coefs")
		guessing_wave[0] = tmp[0]
	else
		//baseline linear
		wave tmp = $("Baseline Coefs")
		guessing_wave[0] = tmp[0]
		guessing_wave[1] = tmp[1]
	endif
	for(i=0;i<num_of_gaussian;i+=1)
		if(hold_wave[0] == 1)
			setdatafolder backup_wd
			wave tmp = $("Peak " + num2str(i) + " Coefs")		//reference to first data (Freq)
			guessing_wave[2+i*3] = tmp[0]
			setdatafolder working_dfr	
		else
			wave tmp = $("Peak " + num2str(i) + " Coefs")
			guessing_wave[2+i*3] = tmp[0]						//reference to previous data (Freq)
		endif
		if(hold_wave[1+i] == 1)					// for Height
			setdatafolder backup_wd
			wave tmp = $("Peak " + num2str(i) + " Coefs")
			guessing_wave[2+i*3+1] = tmp[1]
			setdatafolder working_dfr
		else
			wave tmp = $("Peak " + num2str(i) + " Coefs")
			guessing_wave[2+i*3+1] = tmp[1]
		endif
		if(hold_wave[11+i] == 1)				//for width
			setdatafolder backup_wd
			wave tmp = $("Peak " + num2str(i) + " Coefs")
			guessing_wave[2+i*3+2] = tmp[2]
			setdatafolder working_dfr
		else
			wave tmp = $("Peak " + num2str(i) + " Coefs")
			guessing_wave[2+i*3+2] = tmp[2]
		endif
		
			wpi[i][0] = guessing_wave[2+i*3]		//Freq (Location)
			wpi[i][1] = guessing_wave[2+i*3+1]	//width
			wpi[i][2] = guessing_wave[2+i*3+2]    //Height
			wpi[i][3] = 0
			wpi[i][4] = 0
	endfor
	setdataFolder current_working_dfr
	variable /G MPF2_UserCursors = 1
	CreateCoefWavesFromAutoPeakInfo(currentsetnumber,wpi,"Gauss")		//wanna change Gauss to Lorenztian ? 
	MPF2_PutAutoPeakResultIntoList(currentsetnumber,wpi,1)
	WMHL_OpenAContainer(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Baseline")
	for(i=0;i<num_of_gaussian;i+=1)
		WMHL_OpenAContainer(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Peak " + num2str(i))
	endfor
	
	
	// --------- Lorentzian に変えたい場合 ---------------(ただし, Freq, width,などの値はこれらの変更後に入れた方が良い)
	
	
	if(cmpstr(N_AF_Fitting_Func_Name,"Lorentzian") == 0)
		i = 1
		do
			if(strlen(WMHL_GetItemForRowNumber(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", i)) == 0)
				break
			endif
			if(WMHL_RowIsContainer(graphStr+"#MultiPeak2Panel#P1","MPF2_PeakList",i))
				string prevFunc =  WMHL_GetExtraColumnData(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0, i)
				prevFunc = RemoveEnding(prevFunc, MENU_ARROW_STRING)
				string input_func = "Lorentzian"		
				if(cmpstr(prevFunc,input_func))
					WMHL_ExtraColumnData(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0, i,input_func +MENU_ARROW_STRING, 0)
					MPF2_CoefWaveForListRow(currentsetNumber, i, input_func)
				 	String theItem=""
					if (WMHL_RowIsOpen(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", i))
						theItem = WMHL_GetItemForRowNumber(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", i)
						WMHL_CloseAContainer(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", theItem)	// info for this row has changed. Opening the row re-evaluates the info
						if (strlen(theItem)) // if this was set, then the container was open earlier
							WMHL_OpenAContainer(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", theItem)
						endif
					endif
				endif
			endif
			i+=1
		while(1)
	endif
	if(cmpstr(N_AF_Baseline_Func_Name,"Linear") == 0)
		string  theItem2
		WMHL_ExtraColumnData(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0, 0,"Linear" +MENU_ARROW_STRING, 0)
		MPF2_InfoForBaseline(currentSetNumber,"Linear")
		if(WMHL_RowIsOpen(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0))
			theItem2 = WMHL_GetItemForRowNumber(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0)
			WMHL_CloseAContainer(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", theItem2)
			WMHL_OpenAContainer(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", theItem2)
		endif
	endif
	// --------------------------------------------------
	
	if(cmpstr("Constant",N_AF_Baseline_Func_Name) == 0)
		//baseline const
		WMHL_ExtraColumnData(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0, 1, num2str(guessing_wave[0]), 1)		//because col = input + 2
	else
		//baseline linear
		WMHL_ExtraColumnData(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0, 1, num2str(guessing_wave[0]), 1)
		WMHL_ExtraColumnData(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0, 2, num2str(guessing_wave[1]), 1)
	endif
	for(i=0;i<num_of_gaussian;i+=1)
		variable row_num = WMHL_GetRowNumberForItem(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Peak " + num2str(i))
		if(hold_wave[0] == 1)
			WMHL_ExtraColumnData(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", 1, row_num+1, "Hold", 0,setSelWaveValue=48)
		endif
		if(hold_wave[11+i] == 1)
			WMHL_ExtraColumnData(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", 1, row_num+2, "Hold", 0,setSelWaveValue=48)
		endif
		if(hold_wave[1+i] == 1)
			WMHL_ExtraColumnData(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", 1, row_num+3, "Hold", 0,setSelWaveValue=48)
		endif
	endfor
	HO_MPF2_DoFitButtonProc(currentSetNumber,graphStr)
	
	HO_MPF2_PeakResultsButtonProc(graphStr,currentSetNumber)		//show result
	
	string output_result_name = "Integral_"
	setDataFolder Get_New_AutoFit_DataFolder()
	for(i=0;i<num_of_gaussian;i+=1)
		wave main_w = $(output_result_name + num2str(i))
		wave error_w = $(output_result_name + "Error_" + num2str(i))
		make /O /N=(numpnts(main_w)+1) $(output_result_name + num2str(i))
		make /O /N=(numpnts(error_w)+1) $(output_result_name + "Error_" + num2str(i))
		if(cmpstr(N_AF_Fitting_Func_Name,"Gauss") == 0 )
			//gauss
			wave /T /SDFR=current_working_dfr result_table_for_gauss = $("MPF2_ResultsListWave")
			main_w[numpnts(main_w)-1] = str2num(result_table_for_gauss[i][6])
			error_w[numpnts(error_w)-1] = str2num(replacestring("+/-",result_table_for_gauss[i][7],""))
		else
			//Lorentzian
			wave /T /SDFR=current_working_dfr result_table_for_Loren = $("MPF2_ResultsListWave")
			main_w[numpnts(main_w)-1] = str2num(result_table_for_Loren[i][6])
			error_w[numpnts(error_w)-1] = str2num(replacestring("+/-",result_table_for_Loren[i][7],""))
		endif
	endfor
	
	setdataFolder current_working_dfr
	return 0
End



Function TEST_Hoitakuuuuuuuuuuuuuuu()
	string  x_wavename = "root:Freq"
	string y_wavename = "root:baseline_csv3"
	variable position_of_Acursor = 31
	variable position_of_Bcursor = 58
	variable num_of_gaussian = 1
	wave x_w = $(x_wavename)
	wave y_w = $(y_wavename)
	variable RangeBegin,RangeEnd,RangeReversed
	DFREF saved_dfr = GetDataFolderDFR()
	
	
	MPF2_StartNewMPFit(0,"New Graph",y_wavename,x_wavename,1,0)
	nvar currentSetNumber = root:Packages:MultiPeakFit2:currentSetNumber		//getting  "set number"
	string graphStr = WinName(0,1)
	DFREF working_dfr = $(MPF2_FolderPathFromSetNumber(currentSetNumber))
	//showinfo
	//cursor
	//pcsr()
	string Just_y_wavename = TraceNameList(graphStr,";",1+4)
	Just_y_wavename = Replacestring(";",Just_y_wavename,"")
	RangeBegin = position_of_Acursor 
	RangeEnd = position_of_Bcursor
	//print(Just_y_wavename)
	cursor /A=1 /W=$(graphStr) A $(Just_y_wavename) position_of_Acursor
	cursor /A=1 /W=$(graphStr) B $(Just_y_wavename) position_of_Bcursor
	setDataFolder working_dfr
	MPF2_SetDataPointRange(graphStr,y_w,x_w,position_of_Acursor,position_of_Bcursor,position_of_Acursor)
	Variable/G XPointRangeBegin = RangeBegin
	Variable/G XPointRangeEnd = RangeEnd
	Variable/G XPointRangeReversed = RangeReversed
	
	//Variable setNumber = GetSetNumberFromWinName("EditOrAddPeaksGraph")
	string DFpath = MPF2_FolderPathFromSetNumber(currentSetNumber)
	SVAR YWvName = $(DFpath+":YWvName")
	SVAR XWvName = $(DFpath+":XWvName")
	Wave yw = $YWvName
	Wave/Z xw = $XWvName
	
	string gname = graphStr
	//SVAR gname = graphStr
	//make /O /N=(num_of_gaussian,4) $("Editwpi") /wave=Editwpi_w
	Make /D/N=(num_of_gaussian,5)/O $("W_AutoPeakInfo") /wave=wpi
	//Redimension /N=(DimSize(wpi,0)+1,5) $("W_AutoPeakInfo")
	wpi[0][0] = 4136.6 //Editwpi[i][2]
	wpi[0][1] = 0.15	//Editwpi[i][3]
	wpi[0][2] = 0.019	//Editwpi[i][1]
	wpi[0][3] = 0 			//Editwpi[i][3]/2
	wpi[0][4] = 0			//Editwpi[i][3]/2
	
	variable /G MPF2_UserCursors = 1
	
	//string /G SavedFunctionTypes 
	//SavedFunctionTypes = "Constant;Guass;"
	
	//MPF2_RefreshHoldStrings(gname+"#MultiPeak2Panel") 
	//MPF2_AddPeaksToGraph(currentSetNumber,wpi,0,0,0)
	string peakTypeName = "Gauss"
	//print(graphStr)
	variable BaselineRow = WMHL_GetRowNumberForItem(graphStr+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Baseline")
	string baselinestr = WMHL_GetExtraColumnData(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0, BaselineRow)
	print(baselinestr)
	CreateCoefWavesFromAutoPeakInfo(currentsetnumber,wpi,"Gauss")
	MPF2_PutAutoPeakResultIntoList(currentsetnumber,wpi,1)
	//variable baselinesOpen = WMHL_RowIsOpen(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", BaselineRow+1)
	//print(baselinesOpen)
	//WMHL_AddObject(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "", "Peak "+num2str(45), 1)
	svar FuncListString
	//HO_MPF2_DoFitButtonProc(currentSetNumber,gname)
	WMHL_OpenAContainer(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Baseline")
	WMHL_OpenAContainer(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Peak 0")
	//HO_MPF2_DoFitButtonProc(currentSetNumber,gname)
	variable row_num = 4
	WMHL_ExtraColumnData(gname+"#MultiPeak2Panel#P1", "MPF2_PeakList", 1, row_num, "Hold", 0,setSelWaveValue=48)
	HO_MPF2_DoFitButtonProc(currentSetNumber,gname)
End


Function NAF_PopupWaveSelector_X(event,wavepath,windowName,ctrlName)
	variable event
	string wavepath
	string windowName
	string ctrlName
	DFREF saved_dfr = getDataFolderDFR()
	setdatafolder Get_New_AutoFit_DataFolder()
	SVAR N_AF_X_wave_name
	N_AF_X_wave_name = wavepath
	setdatafolder saved_dfr
End

Function NAF_PopupWaveSelector_Start(event,wavepath,windowName,ctrlName)
	variable event
	string wavepath
	string windowName
	string ctrlName
	DFREF saved_dfr = getDataFolderDFR()
	setdatafolder Get_New_AutoFit_DataFolder()
	SVAR N_AF_Start_wave_name
	N_AF_Start_wave_name = wavepath
	setdatafolder saved_dfr
End

Function NAF_PopupWaveSelector_End(event,wavepath,windowName,ctrlName)
	variable event
	string wavepath
	string windowName
	string ctrlName
	DFREF saved_dfr = getDataFolderDFR()
	setdatafolder Get_New_AutoFit_DataFolder()
	SVAR N_AF_End_wave_name
	N_AF_End_wave_name = wavepath
	setdatafolder saved_dfr
End


Function New_AutoFit_Function()
	core_new_autofit_function__()
End

static Function /DF Get_New_AutoFit_DataFolder()
	DFREF dfr = root:New_autoFit
	if(DataFolderRefStatus(dfr) != 1)
		NewDataFolder /O root:New_autoFit
		dfr = root:New_autoFit
	endif
	return dfr
End

Function NAF_PopupMenuAction_for_autofit(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum	// which item is currently selected (1-based)
	String popStr		// contents of current popup item as string
	DFREF saved_DFR = getDataFolderDFR()
	setdatafolder Get_New_AutoFit_DataFolder()
	svar N_AF_Fitting_Func_Name
	N_AF_Fitting_Func_Name = popStr
	setdatafolder saved_DFR
	print(popStr)
End

Function NAF_PopupMenuAction_for_base(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum	// which item is currently selected (1-based)
	String popStr		// contents of current popup item as string
	DFREF saved_DFR = getDataFolderDFR()
	setdatafolder Get_New_AutoFit_DataFolder()
	svar N_AF_Baseline_Func_Name
	N_AF_Baseline_Func_Name = popStr
	setdatafolder saved_DFR
	print(popStr)
End

// main panel
static Function core_new_autofit_function__()
	string panel_name = "New_AutoFit_Panel"
	if(wintype(panel_name) != 0)
		DoWindow /F $(panel_name)
	else
		NewPanel /W=(50,50,240,700) /K=1 as panel_name
		RenameWindow $S_name, $(panel_name)
		DFREF saved_dfr = getdatafolderDFR()
		setDatafolder Get_New_AutoFit_DataFolder()
		string /G N_AF_Start_wave_name=""
		string /G N_AF_End_wave_name=""
		string /G N_AF_X_wave_name=""
		variable /G N_AF_num_of_peaks = 0
		string /G N_AF_Fitting_Func_Name="Gauss"
		string /G N_AF_Baseline_Func_Name="Constant"
		setDataFolder saved_dfr
		
		GroupBox index_Group, pos={10,5},size={170,70},title="select y_waves(UP:Start, Down:End)"
		//SetVariable setStartwave pos={20,25},size={150,20},value=AF_Start_wave_name,title="s_wave :"
		Button NAF_Popup_Start, pos={20,25},size={150,20}
		MakeButtonIntoWSPopupButton(panel_name,"NAF_Popup_Start","NAF_PopupWaveSelector_Start",options=PopupWS_OptionFloat)
		//SetVariable setEndwave pos={20,49},size={150,20},value=AF_End_wave_name,title="e_wave :"
		Button NAF_Popup_End, pos={20,49},size={150,20}
		MakeButtonIntoWSPopupButton(panel_name,"NAF_Popup_End","NAF_PopupWaveSelector_End",options=PopupWS_OptionFloat)
		//string tmp_str = get_all_data_folder()
		GroupBox index2_Group,pos={10,70},size={170,40},title="select an x-wave"
		Button NAF_Popup_X_wave, pos={20,88},size={150,20}
		MakeButtonIntoWSPopupButton(panel_name,"NAF_Popup_X_wave","NAF_PopupWaveSelector_X",options=PopupWS_OptionFloat)
		checkbox NS_autoFit_checkbox pos={10,110},title="fix on / off (for Freq)"
		//checkbox NS_autoFit_checkbox_area pos={10,102},title = "fix on / off (for Area)"
		SetVariable setwave_name pos={10,130},size={150,20},value=N_AF_num_of_peaks,format="%d",title="num of peaks :"
		string tmp_str="Gauss;Lorentzian;"
		string tmp_folder = "\"" + tmp_str + "\""
		string baseline_str = "Constant;Linear;"
		string list_baseline_str = "\"" + baseline_str + "\""
		PopupMenu NF_select_Fitting_Func,mode=1,Value=#tmp_folder,title="Func",pos={10,150},size={17,15},proc=NAF_PopupMenuAction_for_autofit
		PopupMenu NF_select_baseline_Func,mode=1,Value=#list_baseline_str,pos={100,150},size={17,15},proc=NAF_PopupMenuAction_for_base
		Button draw_button_NS_AutoFit, pos={10,176},size={170,20},proc=NAF_DoFIT_Proc,title="Do Fit", fSize=10
		
		
		variable space = 20
		Groupbox option_group1,pos={10,200},size={170,220},title="select parameters to hold (Height)"
		checkbox NS_autoFit_checkbox_area_1 pos={10,200+space},title="0   gauss for Height"
		checkbox NS_autoFit_checkbox_area_2 pos={10,200+2*space},title="1   gauss for Height"
		checkbox NS_autoFit_checkbox_area_3 pos={10,200+3*space},title="2   gauss for Height"
		checkbox NS_autoFit_checkbox_area_4 pos={10,200+4*space},title="3   gauss for Height"
		checkbox NS_autoFit_checkbox_area_5 pos={10,200+5*space},title="4   gauss for Height"
		checkbox NS_autoFit_checkbox_area_6 pos={10,200+6*space},title="5   gauss for Height"
		checkbox NS_autoFit_checkbox_area_7 pos={10,200+7*space},title="6   gauss for Height"
		checkbox NS_autoFit_checkbox_area_8 pos={10,200+8*space},title="7   gauss for Height"
		checkbox NS_autoFit_checkbox_area_9 pos={10,200+9*space},title="8   gauss for Height"
		checkbox NS_autoFit_checkbox_area_10 pos={10,200+10*space},title="9   gauss for Height"
		
		Groupbox option_group2,pos={10,420},size={170,220},title="---------------------- (width)"
		checkbox NS_autoFit_checkbox_width_1 pos={10,420+space},title="0   gauss for width"
		checkbox NS_autoFit_checkbox_width_2 pos={10,420+2*space},title="1   gauss for width"
		checkbox NS_autoFit_checkbox_width_3 pos={10,420+3*space},title="2   gauss for width"
		checkbox NS_autoFit_checkbox_width_4 pos={10,420+4*space},title="3   gauss for width"
		checkbox NS_autoFit_checkbox_width_5 pos={10,420+5*space},title="4   gauss for width"
		checkbox NS_autoFit_checkbox_width_6 pos={10,420+6*space},title="5   gauss for width"
		checkbox NS_autoFit_checkbox_width_7 pos={10,420+7*space},title="6   gauss for width"
		checkbox NS_autoFit_checkbox_width_8 pos={10,420+8*space},title="7   gauss for width"
		checkbox NS_autoFit_checkbox_width_9 pos={10,420+9*space},title="8   gauss for width"
		checkbox NS_autoFit_checkbox_width_10 pos={10,420+10*space},title="9   gauss for width"
	endif
End


// Do Fit Function (Main Procedure)
Function NAF_DoFIT_Proc(ctrlName) :ButtonControl
	string ctrlName
	DFREF saved_DFR = GetDataFolderDFR()
	setDataFolder Get_New_AutoFit_DataFolder()
	Svar N_AF_Start_wave_name
	Svar N_AF_End_wave_name
	Svar N_AF_X_wave_name
	Nvar N_AF_num_of_peaks
	Svar N_AF_Fitting_Func_Name
	SVar N_AF_Baseline_Func_Name
	if(N_AF_num_of_peaks < 1)
		doalert 0, "num of peaks value has to be more than 1 at least."
		setdatafolder saved_DFR
		return -1
	elseif(N_AF_num_of_peaks > 10)
		doalert 0,"This Function cannot work correctly with more than 10 peaks."
		setdatafolder saved_DFR
		return -2
	endif
	if(0 && cmpstr(N_AF_Fitting_Func_Name,"Lorentzian") == 0)
		doalert 0, "Lorentzian has not been supported yet."
		setdatafolder saved_DFR
		return -3
	endif
	if(cmpstr(N_AF_Start_wave_name,"") == 0)
		doalert 0, "Select start_wave"
		setdatafolder saved_DFR
		return -4
	endif
	if(cmpstr(N_AF_End_wave_name,"") == 0)
		doalert 0, "Select end_wave"
		setdatafolder saved_DFR
		return -5
	endif
	if(cmpstr(N_AF_X_wave_name,"") == 0)
		doalert 0,"Select X wave"
		setdatafolder saved_DFR
		return -6
	endif	
	//get data names(start -> end)
	wave /Z start_w = $(N_AF_Start_wave_name)
	wave /Z end_w = $(N_AF_End_wave_name)
	wave /Z X_wave = $(N_AF_X_wave_name)
	if(waveexists(start_w) == 0 || waveexists(end_w) == 0 || waveexists(X_wave) == 0)
		doalert 0, "some waves you selected don't exist."
		setdatafolder saved_DFR
		return -7
	endif
	string waves_start_dfr_str = GetWavesDataFolder(start_w,1)
	string waves_end_dfr_str = GetWavesDataFolder(end_w,1)
	if(cmpstr(waves_start_dfr_str,waves_end_dfr_str) != 0)
		doalert 0, "Y-waves you selected don't exist on same directory."
		setdatafolder saved_DFR
		return -8
	endif
	string get_str = ""
	DFREF waves_dfr = $(waves_start_dfr_str)
	setDataFolder waves_dfr
	string input_str = wavelist("*",";","")
	//print(input_str)
	setDataFolder Get_New_AutoFit_DataFolder()
	//splitstring /E=req_exp  input_str, get_str
	variable i
	variable Start_Flag_Local = 0
	for(i=0;i<itemsinlist(input_str);i+=1)
		string ele = stringfromlist(i,input_str)
		if(cmpstr(nameofwave(start_w),ele) == 0)
			Start_Flag_Local = 1
		endif
		if(Start_Flag_Local)
			get_str += ele + ";"
		endif
		if(cmpstr(nameofwave(end_w),ele) == 0)
			break
		endif
	endfor
	print(get_str)
	if(cmpstr(get_str,"") == 0)
		doalert 0, "wave is Nothing"
		setdatafolder saved_DFR
		return -9
	endif
	//Check DataPoint
	variable data_length = numpnts(start_w)
	if(data_length != numpnts(X_wave))
		doalert 0, "The DataLength is defferent between Y-waves to X-Wave."
		setdatafolder saved_DFR
		return -10
	endif
	for(i=1;i<itemsinlist(get_str);i+=1)
		string ele2 = stringfromlist(i,get_str)
		wave /SDFR=waves_dfr tmp_w = $(ele2)
		if(numpnts(tmp_w) != data_length)
			doalert 0, "The DataLength of waves you selected  is different."
			setdatafolder saved_DFR
			return -11
		endif
	endfor
	setDataFolder saved_DFR
	string Valid_DFR = getDataFolder(1)
	string result_str
	string reg_exp = "root:Packages:MultiPeakFit2:(.*)"
	splitstring /E=reg_exp Valid_DFR, result_str
	if(cmpstr(result_str,"") == 0)
		doalert 0,"Current directory isn't MultipeakFit2: ...... "
		setdatafolder saved_DFR
		return -12
	endif
	//make holdingwave to determine whether checkbox is set
	setDataFolder Get_New_AutoFit_DataFolder()
	make /O /N=(21) $("Holding_Flag") /wave=Hold_Wave	//21は1 is Freq, 10 is Height, 10 is width, 0 is false, 1 is true.
	make /O /N=(32) $("Guessing_Wave") /wave=Gussing_Wave 	//const is 1, linear 2, Gaussで使うパラメータ3 x 10個
	Hold_Wave = 0	//For initialization
	Gussing_Wave = 0
	setDataFolder saved_DFR
	string panel_name = "New_AutoFit_Panel"
	controlinfo /W=$(panel_name) NS_autoFit_checkbox
	Hold_Wave[0] = V_Value
	for(i=1;i<11;i+=1)
		controlinfo /W=$(panel_name) $("NS_autoFit_checkbox_area_" + num2str(i))
		Hold_Wave[i] = V_Value
	endfor
	for(i=11;i<21;i+=1)
		controlinfo /W=$(panel_name) $("NS_autoFit_checkbox_width_" + num2str(i-10))
		Hold_Wave[i] = V_Value
	endfor
	
	DFREF backup_saved_DFR = saved_DFR
	//ここでFor分で一つずつ処理をしていく
	variable error_num 
	variable a_index,b_index
	Nvar XPointRangeBegin
	Nvar XPointRangeEnd
	a_index = XPointRangeBegin
	b_index = XPointRangeEnd
	Svar GraphName
	variable BaselineRow = WMHL_GetRowNumberForItem(GraphName+"#MultiPeak2Panel#P1", "MPF2_PeakList", "Baseline")
	string baselinestr = WMHL_GetExtraColumnData(GraphName+"#MultiPeak2Panel#P1", "MPF2_PeakList", 0, BaselineRow)
	baselinestr = replacestring("\\",baselinestr,";")
	baselinestr = stringfromlist(0,baselinestr)
	if(cmpstr(baselinestr,N_AF_Baseline_Func_Name) != 0)
		doalert 0, "Baseline Function is different from one you selected."
		setdatafolder backup_saved_DFR
		return -13
	endif
	string output_result_name = "Integral_"
	setDataFolder Get_New_AutoFit_DataFolder()
	for(i=0;i<N_AF_num_of_peaks;i+=1)
		make /O /N=(0) $(output_result_name + num2str(i))
		make /O /N=(0) $(output_result_name + "Error_" + num2str(i))
	endfor
	setDataFolder saved_DFR
	// --------- Gauss も判定
	for(i=0;i<itemsinlist(get_str);i+=1)
		string ele3 = stringfromlist(i,get_str)
		error_num = NAF_Element_handler(backup_saved_DFR,saved_DFR,Get_New_AutoFit_DataFolder(),waves_dfr,ele3,X_wave,Hold_Wave,Gussing_Wave,N_AF_num_of_peaks,a_index,b_index)	//FUNCTION
		if(error_num != 0)
			doalert 0 ,"Cannot read wave. This Function terminated  accidentaly."
			setdatafolder backup_saved_DFR
			return -14 
		endif
		saved_DFR = getDataFolderDFR()
	endfor
	
	setdataFolder backup_saved_DFR
End

#endif

//what is MPF2_SaveFunctionTypes ?

