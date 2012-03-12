class PageController < ApplicationController
  def home
  end
  
  def userguide
	send_file "app/assets/Doc/PFSP_1.2_Userguide.pdf"
  end

  def PFSPwin
	send_file "app/assets/Doc/PFSP_1.2_Multi_Core.rar"
  end
  
  def methdns
	send_file "app/assets/Doc/perlcode.tgz"
  end
end
