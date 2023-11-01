from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
from argparse import ArgumentParser
import os 

def parse_args():
	parser = ArgumentParser(description = "Porcess some arguments")
	parser.add_argument('-s','--species', required = True, help = "scientific organism/species name")
	return parser.parse_args()

def main():
	args = parse_args()
	name = args.species

	current_path = os.getcwd()
	gecko_path = os.path.join(current_path, '/geckodriver')

	# get download destination path
	options = Options()
	driver = webdriver.Firefox(executable_path=gecko_path, options=options)
	driver.get("https://www.ncbi.nlm.nih.gov/sra?term=%28%22"+name+"%22%5BOrganism%5D%20OR%20"+name+"%5BAll%20Fields%5D%29%20AND%20%28cluster_public%5Bprop%5D%20AND%20%22biomol%20dna%22%5BProperties%5D%20AND%20%22strategy%20wgs%22%5BProperties%5D%20AND%20%22library%20layout%20paired%22%5BProperties%5D%20AND%20%22platform%20illumina%22%5BProperties%5D%20AND%20%22strategy%20wgs%22%5BProperties%5D%20OR%20%22strategy%20wga%22%5BProperties%5D%20OR%20%22strategy%20wcs%22%5BProperties%5D%20OR%20%22strategy%20clone%22%5BProperties%5D%20OR%20%22strategy%20finishing%22%5BProperties%5D%20OR%20%22strategy%20validation%22%5BProperties%5D%20AND%20%22filetype%20fastq%22%5BProperties%5D%29&cmd=DetailsSearch")
	driver.implicitly_wait(1)
	buttons = driver.find_element(By.ID, "sendto")
	buttons.click()
	driver.find_element("xpath","//*[@id='dest_File']").click()
	driver.implicitly_wait(1)
	elem = driver.find_element("xpath","//*[@id='file_format']")
	elem.send_keys("runinfo")
	driver.implicitly_wait(1)
	driver.find_element(By.CLASS_NAME, "button_apply.file.ncbipopper-close-button").click()
	#driver.quit()

	# rename and move file to destination path 

if __name__ == '__main__':
	main()