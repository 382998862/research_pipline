# -*- coding: utf-8 -*- 
'''
This script was used to convert xml config file to html web page.

Author:
       huangls
Version:
        1.0;2015-8-17
'''

import sys, os, argparse, glob, os.path,time

reload(sys)
sys.setdefaultencoding('utf8')
from pyh import *
import numpy as np
import math
import re
import xml.dom.minidom 
from xml.dom.minidom import  parseString, getDOMImplementation
from collections import OrderedDict
import urllib
import shutil
import codecs
import uuid
from PIL import Image
Bin=os.path.split(os.path.realpath(__file__))[0]
###global variable
###xugl 20161021
####for ref mRNA
src_path="BMK_9_html/BMK_2_src"
###for others
#src_path="src"

class parseXML(object): 
    def __init__(self, file): 
        self.f = file
    def __iter__(self): 
        return self 
    def parse(self):
        file=open(self.f,'r')
        
    def next(self): 
        r = [self.s,self.e] 
        self.s, self.e = self.bin+self.s, self.e + self.bin 
        return r 
def cp_template(o,s=None):
    if os.path.exists(o+'/'+src_path):
        sys.stdout.write("output dir has \"src\" dir: %s ,that will be removed!!\n" %(o+'/'+src_path))
        shutil.rmtree(o+'/'+src_path)
    if s:
        shutil.copytree(s,o+'/'+src_path)
    else:
        shutil.copytree(Bin+'/src',o+'/'+src_path)
    
def convert_pic(name,path,out,workDir,wMax=1000.):
    try:
        img = Image.open(path)
        w,h=img.size
        if w >wMax:
            r=wMax/w*100
            sys.stdout.write("start: convert -resize %s%%x%s%% %s %s \n" %(r,r,path,out))
            os.system("convert -resize %s%%x%s%% %s %s " %(r,r,path,out))
        else:
            shutil.copyfile(path,out)
        #print dout+'/src/images/'+pic_name
    except IOError:
        sys.stderr.write("file %s not find ,try %s .....\n"%(path,workDir+'/'+path))
        try:
            img = Image.open(workDir+'/'+path)
            w,h=img.size
            if w >wMax:
                r=wMax/w*100
                sys.stdout.write("start: convert -resize %s%%x%s%% %s %s \n" %(r,r,workDir+'/'+path,out))
                os.system("convert -resize %s%%x%s%% %s %s " %(r,r,workDir+'/'+path,out))
            else:
                shutil.copyfile(workDir+'/'+path,out)
        except Exception,e:
            sys.stderr.write("Error:the path of \"%s\" file %s is blank,please check\n" % (name, path))
            exit(1)
def addTable(tmp,isfirst):
    result=""
    if  isfirst:
        for i in tmp:
            result+="<th>%s</th>"%i
    else:
        for i in tmp:
            
            result+="<td>%s</td>"%i
    return result
    
def timenow():
    """return current time as a string
    """
    return time.strftime('[%d/%m/%Y %H:%M:%S]', time.localtime(time.time()))


parser = argparse.ArgumentParser(description='This script was used to convert xml config file to html web page. \nNote: if any file has Chinese characters must save as UTF-8 format or an error about DecodeError will throw out.')
parser.add_argument('-i', '--input', dest='input', required=True, help='input main xml  file')
parser.add_argument('-t', '--temp', dest='temp', required=False, help='input template src  dir, default:%s'% Bin+"/src")
parser.add_argument('-l', '--tableLineNum', dest='tableLineNum', required=False,type=int,default=20, help='Input max table line numbers to show, default:20')
parser.add_argument('-w', '--width', dest='width', required=False,type=float,default=1500., help='Input max width of pic to show,if a pic\'s width larger than this number , it will be Compressed to be this width, default:1500')
parser.add_argument('-o','--outDir',dest='outDir',required=False,default=os.getcwd(),help='specify the output file dir,default is current dir')
parser.add_argument('-n','--name',dest='name',required=False,default='index',help='specify the output file prefix,default is "index"')
args = parser.parse_args()
###################################################################
if not os.path.exists(args.outDir): os.mkdir(args.outDir)
dout=os.path.abspath(args.outDir)
cp_template(dout,args.temp)

args.input=os.path.abspath(args.input)
workDir=os.path.dirname(args.input)
lineNum=args.tableLineNum
picTypes=["img-width-max", "img-width-normal" , "img-width-min"]
sys.stdout.write("%s:Start Convert xml to html ...\n"%timenow())

##########################################################################
dom = xml.dom.minidom.parse(args.input)
root = dom.documentElement

################################get info#################################################################
report_name = root.getElementsByTagName('report_name')[0].getAttribute('value')
report_version = root.getElementsByTagName('report_version')[0].getAttribute('value')
report_code = root.getElementsByTagName('report_code')[0].getAttribute('value')
report_user = root.getElementsByTagName('report_user')[0].getAttribute('value')
report_user_addr = root.getElementsByTagName('report_user_addr')[0].getAttribute('value')
report_time = root.getElementsByTagName('report_time')[0].getAttribute('value')
report_abstract = root.getElementsByTagName('report_abstract')[0].getAttribute('value')
 
#####################################start page####################################################

page = PyH(report_name)
page<< meta(charset="UTF-8")
page.addAN(ann="<!--[if lt IE 9]>\n<script src=\"%s/js/html5shiv.min.js\"  %(src_path)></script><script src=\"%s/js/respond.min.js\" %(src_path)></script>\n<![endif]-->")

page<< meta(author="huangls@biomarker.com.cn, designed by zhengj@biomarker.com.cn")
page<< meta( http_equiv="X-UA-Compatible" ,content="IE=edge")
page<< meta(name="viewport" ,content="width=device-width, initial-scale=1")
page.addCSS(src_path+'/css/bootstrap.min.css',
            src_path+'/css/index.css',
            src_path+'/js/fancyBox/jquery.fancybox.css',
            src_path+'/css/nav.css',
            src_path+'/css/raxus.css',
            
            )
page.addJS(src_path+'/js/jquery-1.11.3.min.js',
           src_path+'/js/nav.js',
           src_path+'/js/raxus-slider.min.js',
           src_path+'/js/fancyBox/jquery.fancybox.pack.js',
           src_path+'/js/fancyBox/jquery.mousewheel-3.0.6.pack.js',
           src_path+'/js/bootstrap.min.js',
           src_path+'/js/ready.js',
           src_path+'/js/scrolltop.js'
           )
page << title(report_name)


container_shadow=page << div(cl='container shadow') 

head= container_shadow << header()
#header << h2(report_name,cl="text-center", id="title")
head << img( src=src_path+"/images/logo.jpg", cl="pull-right")

rowDiv =container_shadow<< div(cl='row')
main=rowDiv<<div(cl="col-md-9", role="main")
headtitle= main << header()
headtitle << h2(report_name,cl="text-center", id="title")

navMain=rowDiv<<div(cl="col-md-3", role="complementary")

nav= navMain << nav( cl="bs-docs-sidebar hidden-print hidden-xs hidden-sm affix-top")
#navigate
ul1=nav<<ul( cl="nav bs-docs-sidenav")


#parase XML config
dom = xml.dom.minidom.parse(args.input)
root = dom.documentElement
num=0
firstH1=True
if root.nodeName=='report':
    for i in root.childNodes:
        if i.nodeType==1:
            num+=1
            if i.nodeName=='report_version':
                pass
            elif i.nodeName=='report_name':
                pass
            elif i.nodeName=='report_code':
                pass
            elif i.nodeName=='report_user':
                pass
            elif i.nodeName=='report_user_addr':
                pass
            elif i.nodeName=='report_time':
                pass
            elif i.nodeName=='report_abstract':
                
                abstract=i.getAttribute('value').encode("UTF-8")
                if not re.match('^\s*$', abstract):
                
                    sect= main << section()
                    
                    sect<<h3("摘要",id='a%s'%num)
                    sect<<p(abstract,cl='text-center')
                    li_1= ul1 <<li() 
                    li_1<<a("摘要",href='#a%s'%num)
                else:
                    sys.stdout.write("Value of report_abstract is empty，report_abstract will not view on web report\n")
                
            elif i.nodeName=='h1':
                name=i.getAttribute('name').encode("UTF-8")
                type=i.getAttribute('type').encode("UTF-8")
                desc=i.getAttribute('desc').encode("UTF-8")
                #nav
                li_1= ul1 <<li() 
                li_1<<a(name,href='#a%s'%num)
                ul2=li_1<<ul(cl='nav')
                #mainpage
                sect= main << section()
                sect<<h3(name,id='a%s'%num)
#                 if firstH1:
#                     try:
#                         abart=sect<<article()
#                         abart<<p(abstract,cl='text-center')
#                     except NameError:
#                         pass
#                     firstH1=False
                #art= sect << article()
                
                #ulNav
                
            elif i.nodeName=='h2':
                name=i.getAttribute('name').encode("UTF-8")
                type=i.getAttribute('type').encode("UTF-8")
                desc=i.getAttribute('desc').encode("UTF-8")
                li_1_2=ul2<<li()
                li_1_2<<a(name,cl="js-active-bg",href='#a%s'%num)
                ul3=li_1_2<<ul(cl='nav')
                #mainpage
                sect<<h4(name,cl="title-h4" ,id='a%s'%num)
                #sect= art << section()
            elif i.nodeName=='h3':
                name=i.getAttribute('name').encode("UTF-8")
                type=i.getAttribute('type').encode("UTF-8")
                desc=i.getAttribute('desc').encode("UTF-8")
                li_1_3=ul3<<li()
                li_1_3<<a(name,cl="js-active-bg",href='#a%s'%num)
                #mainpage
                sect<<h5(name,cl="title-h5" ,id='a%s'%num)
            elif i.nodeName=='h4':
                pass
            elif i.nodeName=='h5':
                pass
            elif i.nodeName=='h6':
                pass
            elif i.nodeName=='p':
                #print 'p',i,
                type=i.getAttribute('type').encode("UTF-8")
                desc=i.getAttribute('desc').encode("UTF-8")
                sect <<p(desc,cl="paragraph")
            elif i.nodeName=='table':
                name=i.getAttribute('name').encode("UTF-8")
                type=i.getAttribute('type').encode("UTF-8")
                desc=i.getAttribute('desc').encode("UTF-8")
                path=i.getAttribute('path').encode("UTF-8")
                try:
                    sect<<p(name,cl="text-center")
                except NameError:
                    sys.stderr.write("Must has one h1 label,please check\n")
                    exit(1)
                tableDiv=sect<<div(cl="table-responsive")
                if not re.match('^\s*$', desc):sect<<p (desc,cl="table-mark small")
                #mytable=tableDiv<<table( cl="table table-hover table-striped")
                mytable=tableDiv<<table( cl="table table-bordered table-hover table-striped")
                
                tablehead=mytable<<thead()
                tablebody=mytable<<tbody()
                firstline=True
                try:
                    #tb=codecs.open(path,'r','utf-8')
                    tb=open(path,'r')
                except IOError:
                    #tb=codecs.open(workDir+'/'+path,'r','utf-8')
                    try:
                        tb=open(workDir+'/'+path,'r')
                    except Exception,err:
                        sys.stderr.write("Error:\"%s\" donot have file,%s" %(name,err))
                        exit(1)
                if type=='full':
                    for k in tb:
                        
                        k=k.strip('\n')
                        if(re.match("^\s*$", k)):continue
                        tmp=re.split('\t', k)
                        if firstline:
                            line=addTable(tmp,firstline)
                            #tablehead<<tr(line,cl='bg-primary')
                            tablehead<<tr(line,cl='bg-info')
                            firstline=False
                        else:
                            line=addTable(tmp,firstline)
                            tablebody<<tr(line)
                else:
                    n=0
                    for k in tb:
                        n+=1
                        k=k.strip('\n')
                        tmp=re.split('\t', k)
                        if firstline:
                            line=addTable(tmp,firstline)
                            tablehead<<tr(line,cl='bg-info')
                            firstline=False
                        else:
                            line=addTable(tmp,firstline)
                            tablebody<<tr(line)
                        if n>lineNum:
                            break
            elif i.nodeName=='file':
                type=i.getAttribute('type').encode("UTF-8")
                desc=i.getAttribute('desc').encode("UTF-8")
                name=i.getAttribute('name').encode("UTF-8")
                path=i.getAttribute('path').encode("UTF-8")
                try:
                    sect<<p("<a href=\"%s\" title=\"click\" class=\"mylink\" target=\"_blank\">%s</a>" %(path,name),cl="paragraph")
                    if not re.match('^\s*$', desc):sect<<p(desc,cl="table-mark small")
                except NameError:
                    sys.stderr.write("Must has one h1 label,please check\n")
                    exit(1)
            elif i.nodeName=='pic':
                type=i.getAttribute('type').encode("UTF-8")
                desc=i.getAttribute('desc').encode("UTF-8")
                name=i.getAttribute('name').encode("UTF-8")
                path=i.getAttribute('path').encode("UTF-8")
                pic_name=os.path.basename(path)
                if not os.path.exists(dout+'/'+src_path+'/images/'+pic_name):
                    new_name=pic_name
                    convert_pic(name,path,dout+'/'+src_path+'/images/'+pic_name,workDir,args.width)
                else:
                    new_name=str(uuid.uuid4())+"."+pic_name.split('.')[-1]
                    sys.stderr.write("same pic name %s,new name will give %s \n"%(pic_name,new_name))
                    convert_pic(name,path,dout+'/'+src_path+'/images/'+new_name,workDir,args.width)
                    
                pic_p=sect <<p(cl="text-center")
                a_pic=pic_p<<a(cl='img-toggle' ,href=src_path+'/images/'+new_name,title=name)
                
                if (type in picTypes ):
                    a_pic<<img(cl=type, src=src_path+'/images/'+new_name, alt=pic_name)
                else:
                    sys.stdout.write('pic view type of \"%s\" about %s not predefined,\"img-width-max\" will be used\n'%(type,path))
                    a_pic<<img(cl="img-width-max", src=src_path+'/images/'+new_name, alt=pic_name)
                    
                sect<< p(name,cl='img-mark text-center small')
                if not re.match('^\s*$', desc):sect<<p(desc,cl="paragraph-mark center-block small img-width-max")
                
            elif i.nodeName=='pic_list':
                type=i.getAttribute('type').encode("UTF-8")
                desc=i.getAttribute('desc').encode("UTF-8")
                name=i.getAttribute('name').encode("UTF-8")
                imgDiv= sect << div(cl="raxus-slider", data_autoplay="3000", data_arrows="show", data_keypress="true", data_thumbnail="bottom")
                imgUl=imgDiv<<ul(cl="slider-relative")
                try:
                    sect<<p(name,cl="img-mark text-center small")
                except NameError:
                    sys.stderr.write("Must has one h1 label,please check\n")
                    exit(1)
                if not re.match('^\s*$', desc):sect<<p(desc,cl="paragraph-mark center-block small img-width-max")
                imglist=i.getElementsByTagName('pic')
                for j in imglist:
                    path=j.getAttribute('path').encode("UTF-8")
                    desc=j.getAttribute('desc').encode("UTF-8")
                    name=j.getAttribute('name').encode("UTF-8")
                    pic_name=os.path.basename(path)
                    if not os.path.exists(dout+'/'+src_path+'/images/'+pic_name):
                        new_name=pic_name
                        convert_pic(name,path,dout+'/'+src_path+'/images/'+pic_name,workDir,args.width)
                    else:
                        new_name=str(uuid.uuid4())+"."+pic_name.split('.')[-1]
                        sys.stderr.write("same pic name %s,new name will give %s \n"%(pic_name,new_name))
                        convert_pic(name,path,dout+'/'+src_path+'/images/'+new_name,workDir,args.width)
#                         try:
#                             shutil.copyfile(path,dout+'/src/images/'+pic_name)
#                         except IOError:
#                             shutil.copyfile(workDir+'/'+path,dout+'/src/images/'+pic_name)
                    imgLi=imgUl<<li(cl="slide")
                    imgA=imgLi<<a(cl="img-toggle", href=src_path+'/images/'+new_name ,title=name)
                    imgA<<img(src=src_path+'/images/'+new_name)
            elif i.nodeName=='file_list':
                type=i.getAttribute('type').encode("UTF-8")
                desc=i.getAttribute('desc').encode("UTF-8")
                name=i.getAttribute('name').encode("UTF-8")
                try:
                    sect<<p(name,cl="paragraph")
                except NameError:
                    sys.stderr.write("Must has one h1 label,please check\n")
                    exit(1)
                    
                filelist=i.getElementsByTagName('file')
                fileol=sect<<ol()
                for j in filelist:
                    desc=j.getAttribute('desc').encode("UTF-8")
                    path=j.getAttribute('path').encode("UTF-8")
                    name=j.getAttribute('name').encode("UTF-8")
                    fileol<<li("<a href=\"%s\" title=\"click\" class=\"mylink\" target=\"_blank\">%s</a>"%(path,name))
                    
            elif i.nodeName=='ref_list':
                type=i.getAttribute('type').encode("UTF-8")
                desc=i.getAttribute('desc').encode("UTF-8")
                name=i.getAttribute('name').encode("UTF-8")
                ul1<<li("<a  href=\"#js-reference\">%s</a>"%name)
                #ul1<<li("<a class=\"js-active-bg\" href=\"#js-reference\">%s</a>"%name)
                refSect=main << section()
                refSect<<h3(name,id="js-reference")
                refol =refSect<<ol()
                
                reflist=i.getElementsByTagName('ref')
                for j in reflist:
                    id=j.getAttribute('id').encode("UTF-8")
                    link=j.getAttribute('link').encode("UTF-8")
                    name=j.getAttribute('name').encode("UTF-8")
                    refol<<li("<a id=\"ref%s\"  href=\"%s\" title=\"click\" target=\"_blank\">%s</a>"%(id,link,name))
                    
            elif i.nodeName=='sup_list':
                type=i.getAttribute('type').encode("UTF-8")
                desc=i.getAttribute('desc').encode("UTF-8")
                name=i.getAttribute('name').encode("UTF-8")
                
                ul1<<li("<a class=\"js-active-bg\" href=\"#js-appendix\">%s</a>"%name)
                
                supSect=main << section()
                supDiv=supSect<<div(id="js-appendix")
                supDiv<<h3(name,id="js-appendix")
                supul =supDiv<<ul(cl="list-unstyled")
                
                suplist=i.getElementsByTagName('sup')
                for j in suplist:
                    id=j.getAttribute('id').encode("UTF-8")
                    link=j.getAttribute('link').encode("UTF-8")
                    name=j.getAttribute('name').encode("UTF-8")
                    supul<<li("<a href=\"#%s\" title=\"click\" target=\"_blank\"><h4 class=\"title-h4\">%s</h4></a>"%(link,name),cl="paragraph")
                
            else:
                print 'attribute \"%s\" not find in predefined xml,please check\n'% i.nodeName
                
else:
    print 'not find <report> node in xml file,please check\n'
    exit(1)


################################################
year=time.strftime('%Y',time.localtime(time.time()))
foot =container_shadow<< footer() 
divFooter=foot<<div( cl="text-center")
divFooter<<p(u'Copyright © 2009-%s 北京百迈客生物科技有限公司版权所有  京ICP备10042835号'.encode("UTF-8") %year) 
divFooter<<p(u'公司地址：北京市顺义区南法信府前街12号顺捷大厦5层'.encode("UTF-8"))
ulfooter=divFooter<<ul(cl="list-inline list-unstyled")
ulfooter<<li("Tel:400-600-3186")
ulfooter<<li("Fax:010-57045001")
ulfooter<<li("Tel:400-600-3186")
ulfooter<<li("E-mail:<a href=\"mailto:tech@biomarker.com.cn\">tech@biomarker.com.cn</a>")
ulfooter<<li(u"微信: BMK-med".encode("UTF-8"))
ulfooter<<li("<a href=\"http://www.biocloud.cn/external/login/toLogin\" target=\"_blank\">百迈客生物云平台</a>")
ulfooter<<li("<a href=\"http://www.biomarker.com.cn/about-us\" target=\"_blank\">关于我们</a>")

#######################################################################
divgoTop=container_shadow<< div(id="goTop")
divgoTopA=divgoTop<<a(cl="backtotop" ,title="Top")
divgoTopA<<img(src=src_path+"/images/goTop.jpg", cl="back-tip")

######################################################################
page.printOut(dout+'/'+args.name+'.html')
sys.stdout.write("%s:End Convert xml to html \nGood luck !\n"%timenow())

