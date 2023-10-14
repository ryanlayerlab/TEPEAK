from Bio import Entrez
from Bio import SeqIO
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
from optparse import OptionParser
import os 


current_path = os.getcwd()
print(current_path)
gecko_path = current_path + "/geckodriver"
print(gecko_path)
parser = OptionParser()

parser.add_option("-s",
	dest="organism",
	help="scientific organism/species name")

(options_args, args) = parser.parse_args()
# get download destination path
# get orgnaism name
if not options_args.organism:
		parser.error('organism name not given')
options = Options()
driver = webdriver.Firefox(executable_path=gecko_path, options=options)
#name = "gallus gallus"
name = options_args.organism
driver.get("https://www.ncbi.nlm.nih.gov/sra?term=%28%22"+name+"%22%5BOrganism%5D%20OR%20"+name+"%5BAll%20Fields%5D%29%20AND%20%28cluster_public%5Bprop%5D%20AND%20%22biomol%20dna%22%5BProperties%5D%20AND%20%22strategy%20wgs%22%5BProperties%5D%20AND%20%22library%20layout%20paired%22%5BProperties%5D%20AND%20%22platform%20illumina%22%5BProperties%5D%20AND%20%22strategy%20wgs%22%5BProperties%5D%20OR%20%22strategy%20wga%22%5BProperties%5D%20OR%20%22strategy%20wcs%22%5BProperties%5D%20OR%20%22strategy%20clone%22%5BProperties%5D%20OR%20%22strategy%20finishing%22%5BProperties%5D%20OR%20%22strategy%20validation%22%5BProperties%5D%20AND%20%22filetype%20fastq%22%5BProperties%5D%29&cmd=DetailsSearch")
driver.implicitly_wait(1)
buttons = driver.find_element(By.ID, "sendto")
buttons.click()
#driver.find_element_by_xpath("//*[@id='dest_File']").click()
driver.find_element("xpath","//*[@id='dest_File']").click()
driver.implicitly_wait(1)
elem = driver.find_element("xpath","//*[@id='file_format']")
elem.send_keys("runinfo")
driver.implicitly_wait(1)
driver.find_element(By.CLASS_NAME, "button_apply.file.ncbipopper-close-button").click()
#driver.quit()

# rename and move file to destination path 