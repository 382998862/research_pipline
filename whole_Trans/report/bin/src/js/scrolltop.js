/**
 * Created by 北京百迈客生物科技有限公司(zhengj) on 2015/8/1.
 */
$(document).ready(function(){
	//1.首先将#goTopBtn隐藏
	$("#goTop").hide();

	//2.给窗口添加事件：滚动时触发
	$(window).scroll(function(){
		//当滚动条的位置处于距顶部0像素以上时，返回顶部按钮出现，否则消失
		if ($(window).scrollTop()>0){
			$("#goTop").fadeIn(800);
			}
		else{
			$("#goTop").fadeOut(800);
			}
	});
	
	//3.当点击返回顶部按钮后，回到页面顶部位置
	$("#goTop").click(function(){
		$('body,html').animate({scrollTop:0},500);
	});
});
