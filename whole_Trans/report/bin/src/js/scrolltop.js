/**
 * Created by ��������������Ƽ����޹�˾(zhengj) on 2015/8/1.
 */
$(document).ready(function(){
	//1.���Ƚ�#goTopBtn����
	$("#goTop").hide();

	//2.����������¼�������ʱ����
	$(window).scroll(function(){
		//����������λ�ô��ھඥ��0��������ʱ�����ض�����ť���֣�������ʧ
		if ($(window).scrollTop()>0){
			$("#goTop").fadeIn(800);
			}
		else{
			$("#goTop").fadeOut(800);
			}
	});
	
	//3.��������ض�����ť�󣬻ص�ҳ�涥��λ��
	$("#goTop").click(function(){
		$('body,html').animate({scrollTop:0},500);
	});
});
