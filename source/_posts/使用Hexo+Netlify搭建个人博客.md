---
title: 使用Hexo+Netlify搭建个人博客
tag: Hexo
categories: 博客
---
#### 个人博客搭建
##### 准备工作
需要安装Git和Node.js  
- [Git](https://git-scm.com/downloads)的安装  

下载windows最新版，Git-2.21.0-64-bit.exe，
在git使用过程中遇到下面的两个错误  
1、git报错：
```
fatal: bad config line 1 in file C:/Users/JIANGXIAOLIANG/.gitconfig  
```
解决方案：
找到提示的目录，然后删掉.gitconfig文件
然后在重新配置用户名和邮箱，输入下面的命令：
```
$ git config --global user.name "用户名"
$ git config --global user.email "邮箱"
```
我在使用的时候没有重新配置用户好像也可以正常使用  
<!--more-->
2、使用hexo 创建目录时
```
C:\Users\Administrator>hexo init F:/blog
INFO  Cloning hexo-starter https://github.com/hexojs/hexo-starter.git
Cloning into 'F:\blog'...
fatal: unable to access 'https://github.com/hexojs/hexo-starter.git/': error setting certificate verify locations:
  CAfile: E:/Using Software/Git/mingw32/ssl/certs/ca-bundle.crt
  CApath: none
WARN  git clone failed. Copying data instead
WARN  Failed to install dependencies. Please run 'npm install' manually!
```
解决方案
```
git config --system http.sslverify false
```



- [Node.js](https://nodejs.org/en/download/)的安装  

在官网下载最新版本的`node-v10.16.0-x64.msi`，像安装一般软件一样安装Node.js,在安装完以后，将`Node`的路径添加到计算的环境变量中，然后验证Node.js的版本
```
C:\Users\Administrator>node --version
v10.16.0
```
##### 安装Hexo
打开git软件下的`git-bash.exe`终端，在终端中使用npm安装Hexo
```
npm install -g hexo-cli
```
安装完成后创建项目文件夹
```
hexo init F:/blog
# 进入刚刚创建的文件夹
cd F:/blog
```
##### 建站
```
#通过npm完成Hexo初始化
npm install
#网站的雏形已经建好了，可以通过hexo服务器来预览成果
hexo server
# 可以前往 http://localhost:4000/ 访问刚刚建立的最新网站
# 新建博客文章
hexo new “我的最新日志”
# 此时已经可以发现在文件夹./source/_posts下面多了一个我的最新日志.md文件
#下一步，生产静态文件
hexo generate
# 如果hexo服务器还在运行中的话，刷新网页，可以看见刚刚创建的博客文章

```

##### 网站发布前的准备工作
有一个细节值得一提，在默认情况下，Hexo将生成的网站内容存储至public文件夹。鉴于我们不需要对该文件夹的内容进行版本控制，我们可将该文件夹添加至.gitignore文档中:
```
echo "/public" >> .gitignore
```
接下来将内容推送到习惯使用的代码托管服务，一般将其托管到GitHub
新建仓库

首先，在GitHub上新建仓库。为了避免出错，在新建仓库时，请不要在创建Initialize this repository with a README前打勾，Add .gitignore和Add a license处请选择None。

鉴于我们的demo基于Hexo和Netlify，在Repository name处填写hexo_netlify来命名仓库。


打开你的电脑终端，切换至你的项目所在的本地文件夹路径：
```
cd /F/blog
# 初始化仓库
git init
# 该命令将创建一个名为.git的子目录，其中包含了你初始化的Git仓库中所需的文件，这些文件是Git仓库的核心。此时，我们仅作了一个初始化的操作，你的项目文件还未被跟踪。  
# 通过git add 来实现对指定文件的跟踪，然后执行git commit提交:
git add .
git commit -m "initial commit"
# 回到之前我们创建GitHub仓库完成的页面，复制远程仓库链接，在终端输入
git remote add origin git@github.com:jingyu9603/hexo_netlify.git
用以下指令推送本地仓库内容至GitHub
git push origin master
```
在上一步的使用中可能出现：
```
$ ssh -T git@github.com
git@github.com: Permission denied (publickey).
```
解决方法
```
# 创建密钥
ssh-keygen -t rsa -b 4096 -C  "11yj3312@gmail.com"
# 当出现“Enter a file in which to save the key”，按回车键就可以了；出现“nter passphrase (empty for no passphrase)”，也是按回车就可以了
# 开启ssh-agent
eval $(ssh-agent -s)
# 添加密钥到ssh中
 ssh-add ~/.ssh/id_rsa
# 将密钥添加到github的ssh中,复制id_rsa.pub中的内容添加到github count setting的SSH and GPG keys 中
$ cat /c/Users/Administrator/.ssh/id_rsa.pub
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAACAQC5piJ696hRcYzvP0OsGdJHXK+OM1fqUMN2cNSlQJ3GbJYnZn2yuJ5+NTy+tXXtieuqmI4m0b7kB1z0kat8i5r/ML9yk8lvIDekZlKBIiSoW8WcgOGr+mhsJJ6tcCwMOXNOcSYA3dI6lQh+ElVbpRyIffr0jhxZVSH792bxcr6FTu8lh0YvHZD0k/N6NcsrZt6yklFJH+o4yB/mxRdkzdAANNI6TuJad+Jght51Djcsd3lY+4tGqRTYSpohYrLofKREvE9+ZxyWv+Nsf+C9I/CTF6h0J14XRnJo7KNx1sXVsVsmBru5BN0rm0M1IuevbJjmJNsWoT6TYICD2WE6tAgigbvjAx5eCZLdqYXmP8J9TqvXI/heQvKU9BhfJETu8ewij8mWA1vrneorNiofQ9ZiIH0OC5sxJQQWL380GeCvzvP4DTJTj+g8GXvmwW/FnPOEVqliE3nKjzg9GghjBzNu1RHSxl5DxBWi+Rkj0vuQYqmbZ1xDEnrI2fD3GLBFu/ew8jZ/lw7Lk2VSA0oDIaFEhWg9wxSh0j8uUqC6RcSF4idpqIwFYCu8Gu1oP0swEkYc/vsIA5P1aPuqeYrUwbKqCJfiCxwd+SJkKsceNJBZm6g3yBjf3hi8DB/pIIbUk0fhwGZYTo7+YItxZQY4haMetXPClzm93QRg1Cx8134a2Q== 11yj3312@gmail.com

# 在SSH and GPG keys 中点击New SSH key，title可以随便命名，将上面的内容复制到其中
```
通过上面的步骤，项目已经推送大GitHub的master分支下面了，接下来我们对hexo进行一些配置：  
在hexo的根目录下面，用vim打开——config.yml文件：
```
# Deployment
## Docs: https://hexo.io/docs/deployment.html
deploy:
  type: git #部署方式
  repository: git@github.com:jingyu9603/yujing_blog.git #关联github仓库
  branch: run-page #部署分支
```
在这里部署的分支用于后面更新博客使用  
在配置好以后，我们再执行下面的操作：
```
hexo clean #清理各种缓存和旧文件
hexo g #生成静态文件
hexo d # 将public目录同步到GitHub

## 在进行部署时可能会遇到下面的错误
ERROR Deployer not found: gi
## 这是由于我们缺少了一个依赖，安装如下:
npm install hexo-deployer-git --save
```
再执行一下`hexo d`操作，我们会发现在我们github项目中多了一个`run-page`的分支

##### 发布网站
- 注册一个Netlify账号，Netlify提供邮箱注册和包括GitHub、GitLab和Bitbucket在内的第三方验证登陆，在这里以GitHub为例
- 注册好以后，进入页面，点击‘New sit from Git’
- 点击GitHub，关联Netlify和github的仓库，首次关联时点击Authorize netlify同意授权后，Netlify可以有权访问你在GitHub上的仓库内容了，授权完成以后，就可以选择之前建立的hexo_netlify仓库
- 选择run-page分支，在Basic build setting 选项中，两个条目清空，最后点击Deploysit。Netlify就可以构建并发布网站内容了
- 构建完成以后，Netlify会在网站发布成功的同时提供给你一个*.netlify.com为后缀以及随机生成的英文名为前缀的子域名。假如你不喜欢Netlify给你的前缀，点击Site settings,在里面可以修改。

#### 博客优化
##### 安装主题
在[网页上](https://www.zhihu.com/question/24422335)寻找自己喜欢的主题，主题的安装过程很简单，以next主题为例子，在官网有安装步骤：
```
mkdir themes/next 
curl -s https://api.github.com/repos/iissnan/hexo-theme-next/releases/latest | grep tarball_url | cut -d '"' -f 4 | wget -i - -O- | tar -zx -C themes/next --strip-components=1
```
##### 启用主题
修改站点的配置文件`_config.yml`
```
# Extensions
## Plugins: https://hexo.io/plugins/
## Themes: https://hexo.io/themes/
theme: next
```
`next`主题下面还分了四种不同的主题：
- Muse - 默认 Scheme，这是 NexT 最初的版本，黑白主调，大量留白
- Mist - Muse 的紧凑版本，整洁有序的单栏外观
- Pisces - 双栏 Scheme，小家碧玉似的清新
- Gemini - 左侧网站信息及目录，块+片段结构布局

在主题配置文件，搜索 scheme 关键字。 你会看到有四行 scheme 的配置，将你需用启用的 scheme 前面注释 # 去除就可以了


##### 使用Valine插件使得Hexo博客具有评论功能

可参见[SmartSi的blog](http://smartsi.club/add-a-comment-function-to-the-next-theme-of-hexo.html)
###### 获取appid和appkey
- 请先登录或注册 [LeanCloud](https://leancloud.cn/dashboard/login.html#/signup), 进入控制台后点击左下角创建应用，选择免费的开发版即可。注意右上角有几个节点，可以就近选择。
- 应用创建好以后，进入刚刚创建的应用，选择左下角的设置>应用Key，然后就能看到你的appid和appkey了：

###### Hexo中的开启和设置

- 在配置文件中修改代码
```
# Valine
valine:
  enable: true
  appid:       ## 之前获得的id
  appkey:      ## 之前获得key
  notify: false #新的留言是否需要通知
  verify: false #留言是否需要验证
  placeholder: 欢迎留言！在这里说出你的想法！
  avater: mm
  guest_info: nick,mail
  pageSize: 10

```

##### 添加搜索功能
- 安装hexo-generator-searchdb 插件
```
npm install hexo-generator-searchdb --save
```
- 打开 站点配置文件 找到Extensions在下面添加
```
# 搜索
search:
  path: search.xml
  field: post
  format: html
  limit: 10000
```
- 打开 主题配置文件 找到Local search，将enable设置为true

##### 添加阅读全文按钮
如果只想要显示文章的部分内容，只需要在文章的合适位置添加：
```
<!--more-->
```

##### 设置网站缩略图标
我的设置方法是这样的：把你的图片（png或jpg格式，不是favicon.ico）放在themes/next/source/images里，然后打开 主题配置文件 找到favicon，将small、medium、apple_touch_icon三个字段的值都设置成/images/图片名.jpg就可以了，其他字段都注释掉。

##### 头像设置
打开 主题配置文件 找到Sidebar Avatar字段
```
# Sidebar Avatar
avatar: /images/header.jpg
```

##### 设置代码高亮
在站点配置文件中，搜索`highlight`
```
highlight:
  enable: true
  line_number: true
  auto_detect: true
  tab_replace:
```
文字自动检测默认不启动，所以改成true使其起作用。  
再到主题配置文件中：
找到`highlight_theme: normal`，注释显示有五种显示主题可用，分别是：
- normal
- night
- night eighties
- night blue
- night bright
选择什么要看个人审美了。